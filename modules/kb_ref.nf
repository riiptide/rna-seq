process KB_REF {
    tag 'Ensembl108-v1'
    label 'light'
    storeDir params.reference_cache
    input:
    path downloads_manifest
    path inputs_manifest
    val url
    val retries
    path code
    output:
    tuple path('index.idx'), path('t2g.txt'), emit: reference
    script:
    """
    python reference.py --url '${url}' --retries ${retries} \
      --downloads-manifest '${downloads_manifest}' --inputs-manifest '${inputs_manifest}'
    """
}
