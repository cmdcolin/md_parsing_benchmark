import fs from 'fs'
import type { Buffer } from 'buffer'
import { Bench } from 'tinybench'
import { ALL } from 'dns'

const mdRegex = /(\d+)|(\^[A-Za-z]+)|([A-Za-z])/g

// get relative reference sequence positions for positions given relative to
// the read sequence
export function* getNextRefPos(cigarOps: string[], positions: number[]) {
  let readPos = 0
  let refPos = 0
  let currPos = 0

  for (let i = 0; i < cigarOps.length && currPos < positions.length; i += 2) {
    const len = +cigarOps[i]!
    const op = cigarOps[i + 1]!
    if (op === 'S' || op === 'I') {
      for (let i = 0; i < len && currPos < positions.length; i++) {
        if (positions[currPos] === readPos + i) {
          currPos++
        }
      }
      readPos += len
    } else if (op === 'D' || op === 'N') {
      refPos += len
    } else if (op === 'M' || op === 'X' || op === '=') {
      for (let i = 0; i < len && currPos < positions.length; i++) {
        if (positions[currPos] === readPos + i) {
          yield {
            readPos: readPos + i,
            refPos: refPos + i,
            idx: currPos,
          }
          currPos++
        }
      }
      readPos += len
      refPos += len
    }
  }
}

export function getNextRefPosNoGen(cigarOps: string[], positions: number[]) {
  let readPos = 0
  let refPos = 0
  let currPos = 0
  const ret = [] as { refPos: number; readPos: number; idx: number }[]
  for (let i = 0; i < cigarOps.length && currPos < positions.length; i += 2) {
    const len = +cigarOps[i]!
    const op = cigarOps[i + 1]!
    if (op === 'S' || op === 'I') {
      for (let i = 0; i < len && currPos < positions.length; i++) {
        if (positions[currPos] === readPos + i) {
          currPos++
        }
      }
      readPos += len
    } else if (op === 'D' || op === 'N') {
      refPos += len
    } else if (op === 'M' || op === 'X' || op === '=') {
      for (let i = 0; i < len && currPos < positions.length; i++) {
        if (positions[currPos] === readPos + i) {
          ret.push({
            readPos: readPos + i,
            refPos: refPos + i,
            idx: currPos,
          })
          currPos++
        }
      }
      readPos += len
      refPos += len
    }
  }
  return ret
}

export function mdToMismatchesYesGen(
  mdstring: string,
  ops: string[],
  seq: string,
  qual?: Buffer | number[] | null,
): any[] {
  // eslint-disable-next-line @typescript-eslint/no-non-null-assertion
  const mdMatches = mdstring.match(mdRegex)!
  const positions = [] as number[]
  const bases = [] as string[]
  let currPos = 0
  for (const match of mdMatches) {
    const token = match
    const c = token[0]
    if (c >= '0' && c <= '9') {
      currPos += +token
    } else if (!token.startsWith('^')) {
      bases.push(token)
      positions.push(currPos)
      currPos += 1
    }
  }
  const res = [] as {
    start: number
    base: string
    length: number
    qual: number | undefined
    type: string
    altbase: string
  }[]
  let i = 0
  for (const { refPos, readPos } of getNextRefPosNoGen(ops, positions)) {
    res.push({
      start: refPos,
      base: seq[readPos]!,
      length: 1,
      qual: qual?.[readPos],
      type: 'mismatch' as const,
      altbase: bases[i],
    })
    i++
  }
  return res
}
export function mdToMismatchesNoGen(
  mdstring: string,
  ops: string[],
  seq: string,
  qual?: Buffer | number[] | null,
): any[] {
  // eslint-disable-next-line @typescript-eslint/no-non-null-assertion
  const mdMatches = mdstring.match(mdRegex)!
  const positions = [] as number[]
  const bases = [] as string[]
  let currPos = 0
  for (const match of mdMatches) {
    const token = match
    const c = token[0]
    if (c >= '0' && c <= '9') {
      currPos += +token
    } else if (!token.startsWith('^')) {
      bases.push(token)
      positions.push(currPos)
      currPos += 1
    }
  }
  const res = [] as {
    start: number
    base: string
    length: number
    qual: number | undefined
    type: string
    altbase: string
  }[]
  let i = 0
  for (const { refPos, readPos } of getNextRefPos(ops, positions)) {
    res.push({
      start: refPos,
      base: seq[readPos]!,
      length: 1,
      qual: qual?.[readPos],
      type: 'mismatch' as const,
      altbase: bases[i],
    })
    i++
  }
  return res
}

/**
 * parse a SAM MD tag to find mismatching bases of the template versus the
 * reference @returns array of mismatches and their positions
 */
export function mdToMismatchesOld(
  mdstring: string,
  ops: string[],
  cigarMismatches: any[],
  seq: string,
  qual?: Buffer | number[] | null,
) {
  let curr: any = { start: 0, base: '', length: 0, type: 'mismatch' }
  let lastCigar = 0
  let lastTemplateOffset = 0
  let lastRefOffset = 0
  let lastSkipPos = 0
  const mismatchRecords: any[] = []
  const skips = cigarMismatches.filter(cigar => cigar.type === 'skip')

  // convert a position on the reference sequence to a position on the template
  // sequence, taking into account hard and soft clipping of reads
  function nextRecord(): void {
    mismatchRecords.push(curr)

    // get a new mismatch record ready
    curr = {
      start: curr.start + curr.length,
      length: 0,
      base: '',
      type: 'mismatch',
    }
  }

  function getTemplateCoordLocal(refCoord: number): number {
    let templateOffset = lastTemplateOffset
    let refOffset = lastRefOffset
    for (
      let i = lastCigar;
      i < ops.length && refOffset <= refCoord;
      i += 2, lastCigar = i
    ) {
      const len = +ops[i]!
      const op = ops[i + 1]!

      if (op === 'S' || op === 'I') {
        templateOffset += len
      } else if (op === 'D' || op === 'P' || op === 'N') {
        refOffset += len
      } else if (op !== 'H') {
        templateOffset += len
        refOffset += len
      }
    }
    lastTemplateOffset = templateOffset
    lastRefOffset = refOffset

    return templateOffset - (refOffset - refCoord)
  }

  // now actually parse the MD string
  const md = mdstring.match(mdRegex) || []
  for (const token of md) {
    const num = +token
    if (!Number.isNaN(num)) {
      curr.start += num
    } else if (token.startsWith('^')) {
      curr.start += token.length - 1
    } else {
      // mismatch
      // eslint-disable-next-line @typescript-eslint/prefer-for-of
      for (let j = 0; j < token.length; j += 1) {
        curr.length = 1

        while (lastSkipPos < skips.length) {
          const mismatch = skips[lastSkipPos]!
          if (curr.start >= mismatch.start) {
            curr.start += mismatch.length
            lastSkipPos++
          } else {
            break
          }
        }
        const s = getTemplateCoordLocal(curr.start)
        curr.base = seq[s] || 'X'
        curr.qual = qual?.[s]
        curr.altbase = token
        nextRecord()
      }
    }
  }
  return mismatchRecords
}

export function cigarToMismatches(
  ops: string[],
  seq?: string,
  ref?: string,
  qual?: Buffer,
) {
  let roffset = 0 // reference offset
  let soffset = 0 // seq offset
  const mismatches: any[] = []
  const hasRefAndSeq = ref && seq
  for (let i = 0; i < ops.length; i += 2) {
    const len = +ops[i]!
    const op = ops[i + 1]!

    if (op === 'M' || op === '=' || op === 'E') {
      if (hasRefAndSeq) {
        for (let j = 0; j < len; j++) {
          if (
            // @ts-ignore in the full yarn build of the repo, this says that
            // object is possibly undefined for some reason, ignored
            seq[soffset + j].toUpperCase() !== ref[roffset + j].toUpperCase()
          ) {
            mismatches.push({
              start: roffset + j,
              type: 'mismatch',
              base: seq[soffset + j]!,
              altbase: ref[roffset + j]!,
              length: 1,
            })
          }
        }
      }
      soffset += len
    }
    if (op === 'I') {
      mismatches.push({
        start: roffset,
        type: 'insertion',
        base: `${len}`,
        length: 0,
      })
      soffset += len
    } else if (op === 'D') {
      mismatches.push({
        start: roffset,
        type: 'deletion',
        base: '*',
        length: len,
      })
    } else if (op === 'N') {
      mismatches.push({
        start: roffset,
        type: 'skip',
        base: 'N',
        length: len,
      })
    } else if (op === 'X') {
      const r = seq?.slice(soffset, soffset + len) || []
      const q = qual?.subarray(soffset, soffset + len) || []

      for (let j = 0; j < len; j++) {
        mismatches.push({
          start: roffset + j,
          type: 'mismatch',
          base: r[j]!,
          qual: q[j]!,
          length: 1,
        })
      }
      soffset += len
    } else if (op === 'H') {
      mismatches.push({
        start: roffset,
        type: 'hardclip',
        base: `H${len}`,
        cliplen: len,
        length: 1,
      })
    } else if (op === 'S') {
      mismatches.push({
        start: roffset,
        type: 'softclip',
        base: `S${len}`,
        cliplen: len,
        length: 1,
      })
      soffset += len
    }

    if (op !== 'I' && op !== 'S' && op !== 'H') {
      roffset += len
    }
  }
  return mismatches
}
const cigarRegex = new RegExp(/([MIDNSHPX=])/)
export function parseCigar(cigar = '') {
  return cigar.split(cigarRegex).slice(0, -1)
}

const bench = new Bench({ name: 'simple benchmark', time: 100 })

const lines = fs
  .readFileSync('out.small.txt', 'utf8')
  .split('\n')
  .map(r => {
    const [seq, md, cigar] = r.split('\t')
    const ops = parseCigar(cigar)
    return {
      seq,
      md,
      ops,
      mismatches: cigarToMismatches(ops),
    }
  })

bench
  .add('generator', () => {
    for (const { seq, md, ops } of lines) {
      mdToMismatchesNoGen(md, ops, seq)
    }
  })
  .add('no generator', () => {
    for (const { seq, md, ops } of lines) {
      mdToMismatchesYesGen(md, ops, seq)
    }
  })
  .add('old generator', () => {
    for (const { seq, md, ops, mismatches } of lines) {
      mdToMismatchesOld(md, ops, mismatches, seq)
    }
  })

await bench.run()

console.log(bench.name)
console.table(bench.table())
