# sourmash_plugin_betterplot

## Installation

```
pip install sourmash_plugin_betterplot
```

## Usage

@CTB: info goes here!

## Examples

```
sourmash scripts plot2 3sketches.cmp -o examples/plot2.3sketches.cmp.png
```
yields:
![examples/plot2.3sketches.cmp.png](examples/plot2.3sketches.cmp.png)

## Support

We suggest filing issues in [the main sourmash issue tracker](https://github.com/dib-lab/sourmash/issues) as that receives more attention!

## Dev docs

`betterplot` is developed at https://github.com/sourmash-bio/sourmash_plugin_betterplot.

### Testing

Run:
```
pytest tests
```

### Generating a release

Bump version number in `pyproject.toml` and push.

Make a new release on github.

Then pull, and:

```
python -m build
```

followed by `twine upload dist/...`.
