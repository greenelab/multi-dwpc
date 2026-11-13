FROM python:3.13

# install uv
RUN pip install uv

# copy in pyproject.toml and install via
COPY pyproject.toml .
RUN --mount=type=cache,target=/root/.cache/ \
    uv pip install --system -r pyproject.toml

WORKDIR /app

# copy in the rest of the code
COPY . .

# this launches the DRF API rather than the streamlit app
WORKDIR /app/api
CMD ["./launch_api.sh"]
