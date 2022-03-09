export const GALEF_PATH = "galef.json";
export async function get(path = GALEF_PATH) {
  const res = await fetch(path);
  if (res.ok) {
    const data = await res.json();
    let idx = 1;
    for (const block of data) {
      for (const q of block.questions) {
        q.idx = idx++;
      }
    }
    return data;
  }
  throw new Error(`error: ${res.status} ${res.statusText}`);
}
