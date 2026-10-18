FROM python:3.11-slim AS builder
COPY --from=ghcr.io/astral-sh/uv:0.13.0 /uv /bin/uv

ENV UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    UV_PYTHON_DOWNLOADS=never

WORKDIR /app

# install dependencies first so they are cached independently of source changes
COPY pyproject.toml uv.lock README.md LICENSE.txt ./
RUN --mount=type=cache,target=/root/.cache/uv uv sync --locked --no-dev --no-install-project

COPY edgehog ./edgehog
RUN --mount=type=cache,target=/root/.cache/uv uv sync --locked --no-dev --no-editable


FROM python:3.11-slim AS runtime
WORKDIR /app

ENV VIRTUAL_ENV=/app/.venv \
    PATH="/app/.venv/bin:$PATH"

COPY --from=builder ${VIRTUAL_ENV} ${VIRTUAL_ENV}

ENTRYPOINT ["edgehog"]
