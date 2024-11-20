FROM python:3.11-slim AS builder


ENV POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_IN_PROJECT=1 \
    POETRY_VIRTUALENVS_CREATE=1 \
    POETRY_CACHE_DIR=/tmp/poetry_cache

RUN pip install poetry
RUN python -m venv /venv

WORKDIR /app

COPY pyproject.toml README.md ./
RUN --mount=type=cache,target=$POETRY_CACHE_DIR poetry lock && poetry install --without dev --no-root


FROM python:3.11-slim AS runtime
WORKDIR /app

ENV VIRTUAL_ENV=/app/.venv \
    PATH="/app/.venv/bin:$PATH"

COPY pyproject.toml README.md ./
COPY --from=builder ${VIRTUAL_ENV} ${VIRTUAL_ENV}
COPY edgehog ./edgehog
RUN pip install poetry && poetry install --without dev && pip uninstall -y poetry

ENTRYPOINT ["edgehog"]
