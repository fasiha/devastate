export const GALEF_PATH = "/galef.json";
export async function get(path = GALEF_PATH) {
  const res = await fetch(path);
  if (res.ok) {
    return res.json();
  }
  throw new Error(`error: ${res.status} ${res.statusText}`);
}
