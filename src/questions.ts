export interface Question {
  question?: string;
  options: [string, string], answer: number;
  comment?: string;
  idx: number;
}

export interface QuestionBlock {
  title: string;
  questions: Question[];
}

export const GALEF_PATH = 'galef.json';

export async function get(path = GALEF_PATH): Promise<QuestionBlock[]> {
  const res = await fetch(path);
  if (res.ok) {
    // TODO validate via io-ts, etc.
    const data: QuestionBlock[] = await res.json();
    let idx = 1;
    for (const block of data) {
      for (const q of block.questions) { q.idx = idx++; }
    }
    return data;
  }
  // TODO report this gracefully
  throw new Error(`error: ${res.status} ${res.statusText}`);
}