# Multi-DWPC frontend

React and TypeScript frontend for the API in `../api`, matching the Streamlit app in `../app.py`.

## Development

From the repo root, this starts the API and the frontend in Docker:

```bash
./run_api.sh
```

Open http://localhost:5173. Changes under `src/` reload automatically. After changing `package.json` or `vite.config.ts`, rebuild by running the script again.

To run the frontend outside Docker instead, start only the API with `docker compose up streamlit-api`, then:

```bash
npm install
npm run dev
```

Vite proxies `/api` to `API_PROXY_TARGET`, which defaults to `http://localhost:8500`.

## Production build

```bash
npm run build
```

This writes static files to `dist/`. The app calls `/api` with relative URLs, so serve it from the same origin as the API.

## API types

`src/api/schema.d.ts` is generated from the API's OpenAPI schema. With the API running, refresh it after API changes:

```bash
npm run gen:api
```
