# POST a form-encoded request to the STRING API

Sends `params` as an `application/x-www-form-urlencoded` POST body (the
same key=value pairs STRING's GET endpoints accept as query parameters,
just relocated to the request body). Used by
[`APOLLO_get_string_ppi`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_get_string_ppi.md)
instead of GET so identifier lists in the hundreds/thousands don't hit
the URL-length ceiling GET requests run into. Values already containing
percent-encoded characters (e.g. the `identifiers` field, joined with
`%0d`) are passed through as-is.

## Usage

``` r
.apollo_string_api_post(path, params, api_base = "https://string-db.org/api")
```

## Arguments

- path:

  Character. API path after `api_base/`, e.g. `"tsv/network"` or
  `"image/network"`.

- params:

  Named list of POST body fields (character/numeric scalars).

- api_base:

  Character. STRING API base URL.

## Value

The raw response list from
[`curl::curl_fetch_memory()`](https://jeroen.r-universe.dev/curl/reference/curl_fetch.html)
(`$status_code`, `$content`, ...).
