import react from "@vitejs/plugin-react";
import { defineConfig } from "vite";

// Requests use relative /api URLs, in development and in production.
const proxy = {
  "/api": {
    target: process.env.API_PROXY_TARGET ?? "http://localhost:8500",
    changeOrigin: true,
  },
};

export default defineConfig({
  plugins: [react()],
  server: { proxy },
  preview: { proxy },
});
