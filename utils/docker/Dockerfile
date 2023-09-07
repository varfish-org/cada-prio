# syntax=docker/dockerfile:1.3

# ---------------------------------------------------------------------------
# Base
# ---------------------------------------------------------------------------

FROM python:3.11-slim AS base

ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONFAULTHANDLER 1

RUN apt-get update && apt-get install -y curl
RUN apt-get clean autoclean && \
    apt-get autoremove --yes && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/


# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------

FROM base AS deps

RUN pip install pipenv
RUN apt-get update && apt-get install -y --no-install-recommends gcc

COPY requirements requirements
RUN virtualenv /.venv && \
    . /.venv/bin/activate && \
    pip install -r requirements/prod.txt

# ---------------------------------------------------------------------------
# Runtime
# ---------------------------------------------------------------------------

FROM base AS runtime

COPY --from=deps /.venv /.venv

ENV PATH="/.venv/bin:$PATH"

COPY cada_prio cada_prio
COPY requirements requirements
COPY setup.py requirements.txt README.md CHANGELOG.md .
RUN pip install .

RUN useradd --create-home cada_prio
WORKDIR /home/cada_prio
USER cada_prio

CMD ["uvicorn", "cada_prio.server:app", "--host", "0.0.0.0", "--port", "8080"]
EXPOSE 8080
