export interface Question {
  question: string;
  answer: boolean;
  comment?: string;
}

export interface QuestionBlock {
  title: string;
  questions: Question[];
}

export const GALEF_PATH = '/galef.json';

export async function get(path = GALEF_PATH): Promise<QuestionBlock[]> {
  const res = await fetch(path);
  if (res.ok) {
    // TODO validate via io-ts, etc.
    return res.json();
  }
  // TODO report this gracefully
  throw new Error(`error: ${res.status} ${res.statusText}`);
}