import { defineConfig } from "astro/config";
import mdx from "@astrojs/mdx";
import tailwind from "@astrojs/tailwind";
import remarkGfm from "remark-gfm";

export default defineConfig({
  site: "https://github.com/Trhova/Transcriptomics-extended",
  integrations: [mdx(), tailwind({ applyBaseStyles: false })],
  markdown: {
    remarkPlugins: [remarkGfm]
  }
});
