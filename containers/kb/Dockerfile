# The upstream image omits two shared libraries required by its bundled tools.
FROM rna-seq-analysis:0.3.0 AS runtime
FROM quay.io/biocontainers/kb-python:0.30.2--pyh106432d_0@sha256:185589a145027d70fcbc9d6ff8fab5fc7651e1af7099b31646973fa2c55fc82a
COPY --from=runtime /usr/lib/x86_64-linux-gnu/libstdc++.so.6* /usr/lib/x86_64-linux-gnu/
COPY --from=runtime /usr/lib/x86_64-linux-gnu/libz.so.1* /usr/lib/x86_64-linux-gnu/
ARG KALLISTO_URL=https://github.com/pachterlab/kallisto/releases/download/v0.51.1/kallisto_linux-v0.51.1.tar.gz
ARG KALLISTO_SHA256=943021dcf2610dd7438604c6b907b64a0765455f7d96cd63de897a9d15d45c51
RUN python -c "import urllib.request,hashlib,tarfile; p='/tmp/kallisto.tar.gz'; urllib.request.urlretrieve('${KALLISTO_URL}',p); assert hashlib.sha256(open(p,'rb').read()).hexdigest() == '${KALLISTO_SHA256}'; t=tarfile.open(p); open('/usr/local/bin/kallisto','wb').write(t.extractfile('kallisto/kallisto').read())" && chmod +x /usr/local/bin/kallisto
RUN kallisto version && /usr/local/lib/python3.12/site-packages/kb_python/bins/linux/bustools/bustools version
