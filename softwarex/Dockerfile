FROM python:3.12-slim

WORKDIR /app
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
COPY examples ./examples

RUN python -m pip install --no-cache-dir .

ENTRYPOINT ["ssi-cert"]
CMD ["verify", "examples/certs/cert_demo.json"]
