FROM node:22-alpine

WORKDIR /app

# install dependencies first, so they stay cached until the lockfile changes
COPY package.json package-lock.json ./
RUN npm ci

# copy in the rest of the code
COPY . .

EXPOSE 5173

# --host makes the dev server reachable from outside the container
CMD ["npm", "run", "dev", "--", "--host", "0.0.0.0"]
