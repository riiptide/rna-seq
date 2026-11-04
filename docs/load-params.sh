# Source from the repository root in Bash or Zsh after activating the environment.
# Nextflow parses its own configuration; Python only converts JSON to quoted shell values.
mkdir -p work
nextflow -C nextflow.config config -o json > work/params.json || return 1
pbmc_exports="$(python - <<'PY'
import json
import shlex
with open('work/params.json') as handle:
    params = json.load(handle)['params']
for key, value in params.items():
    if not key.isidentifier() or not key.isascii():
        raise ValueError(f'Invalid parameter name: {key}')
    if not isinstance(value, (str, int, float, bool)):
        raise ValueError(f'Expected a scalar parameter: {key}')
    print(f'export {key.upper()}={shlex.quote(str(value))}')
PY
)" || return 1
eval "$pbmc_exports"
unset pbmc_exports
