import { defineCollection, z } from "astro:content";

const guides = defineCollection({
  type: "content",
  schema: z.object({
    title: z.string(),
    description: z.string(),
    slug: z.string(),
    tags: z.array(z.string()),
    toc: z.boolean().default(true),
    draft: z.boolean().default(false)
  })
});

export const collections = { guides };
