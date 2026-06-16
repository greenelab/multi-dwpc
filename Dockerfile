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

# run app.py as a streamlit app
CMD ["uv", "run", "app.py", "--host", "0.0.0", "--port", "8501"]
