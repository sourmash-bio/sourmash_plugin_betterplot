# sourmash_plugin_xyz: a template for sourmash plugins

This is a good place to start if you're writing a plugin for
[sourmash (sourmash-bio/sourmash/)](https://github.com/sourmash-bio/sourmash/).

Note: plugins were added in [sourmash v4.7.0](https://github.com/sourmash-bio/sourmash/releases/tag/v4.7.0).

## Instructions

You can use this repo as a template repo to create a new plugin!

See [this set of changes](https://github.com/ctb/sourmash_plugin_template_test1/pull/1) for the minimal diff needed to get a plugin working!

### Building from a template:

First, go to [the GitHub page](https://github.com/sourmash-bio/sourmash_plugin_template) and click "Use this template" to create a new repository.

Clone that repository into your development space.

Then, search for all places where 'xyz' is present, and replace
'xyz' with the name of your plugin.

Next, edit the code in `src/sourmash_plugin_xyz.py` to implement the plugin
(you'll probably want to change the name of that file, too.)

Then run `pip install .` to install and test your plugin! You can also
run `pip install -e .` to install it in editable mode, which is more
convenient for development.

## Examples

[sourmash_plugin_avro](https://github.com/sourmash-bio/sourmash_plugin_avro)
and
[sourmash_plugin_load_urls](https://github.com/sourmash-bio/sourmash_plugin_load_urls)
are two examples you can follow.

## Template docs for new plugin built from this template.

Delete everything from this line on up and put in your new README ;).

# sourmash_plugin_xyz

## Installation

```
pip install sourmash_plugin_xyz
```

## Usage

non-xyz info goes here!

## Support

We suggest filing issues in [the main sourmash issue tracker](https://github.com/dib-lab/sourmash/issues) as that receives more attention!

## Dev docs

`xyz` is developed at https://github.com/sourmash-bio/sourmash_plugin_template.

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
