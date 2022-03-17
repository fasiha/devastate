import npq from './npq.json';

export function get(): Record<number, Record<number, Record<number, number>>> { return npq; }

export function lookup(n: number, p: number, q: number): number {
  const foo = get();
  const ret = foo?.[n]?.[p]?.[q];
  if (typeof ret !== 'number') { throw new Error('unknown n/p/q'); }
  return ret;
}

export function lookups(n: number, p: number): [number, number][] {
  const foo = get();
  const dict = foo?.[n]?.[p];
  if (!dict) { throw new Error('unknown n/p'); }
  return Object.entries(dict).map(([k, v]) => [Number(k), v]);
}