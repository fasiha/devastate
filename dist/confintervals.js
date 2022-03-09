import npq from "./npq.json.proxy.js";
export function get() {
  return npq;
}
export function lookup(n, p, q) {
  const foo = get();
  const ret = foo?.[n]?.[p]?.[q];
  if (typeof ret !== "number") {
    throw new Error("unknown n/p/q");
  }
  return ret;
}
export function lookups(n, p) {
  const foo = get();
  const dict = foo?.[n]?.[p];
  if (!dict) {
    throw new Error("unknown n/p");
  }
  return Object.entries(dict).map(([k, v]) => [Number(k), v]);
}
