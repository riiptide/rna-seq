# The upstream image omits two shared libraries required by its bundled tools.
FROM rna-seq-analysis:0.3.0 AS runtime
FROM quay.io/biocontainers/kb-python:0.30.2--pyh106432d_0@sha256:185589a145027d70fcbc9d6ff8fab5fc7651e1af7099b31646973fa2c55fc82a
COPY --from=runtime /usr/lib/x86_64-linux-gnu/libstdc++.so.6* /usr/lib/x86_64-linux-gnu/
COPY --from=runtime /usr/lib/x86_64-linux-gnu/libz.so.1* /usr/lib/x86_64-linux-gnu/
RUN /usr/local/lib/python3.12/site-packages/kb_python/bins/linux/kallisto/kallisto version && /usr/local/lib/python3.12/site-packages/kb_python/bins/linux/bustools/bustools version
