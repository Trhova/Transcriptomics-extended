/** @type {import('tailwindcss').Config} */
export default {
  content: ["./src/**/*.{astro,html,js,jsx,md,mdx,svelte,ts,tsx,vue}"],
  theme: {
    extend: {
      fontFamily: {
        display: ["'Space Grotesk'", "'IBM Plex Sans'", "sans-serif"],
        sans: ["'IBM Plex Sans'", "'Avenir Next'", "'Segoe UI'", "sans-serif"],
        mono: ["'IBM Plex Mono'", "'SFMono-Regular'", "monospace"]
      },
      colors: {
        ink: {
          950: "#06111a",
          900: "#0c1721",
          800: "#132332",
          700: "#173246",
          300: "#9ab9c8",
          200: "#c9dbe5",
          100: "#edf5f8"
        },
        accent: {
          500: "#7dd3c3",
          400: "#99f6e4",
          300: "#b7ffe6"
        },
        amberish: {
          400: "#f4c680",
          300: "#f8ddb0"
        }
      },
      boxShadow: {
        panel: "0 20px 80px rgba(0, 0, 0, 0.28)"
      }
    }
  },
  plugins: [require("@tailwindcss/typography")]
};
