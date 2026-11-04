#!/usr/bin/env python3
"""Create the old:new folder pairs needed to run the pipeline offline.

Every absolute path in the config that is copied into the design/reference folder by
`hydra-genetics references download` (matched by url in the validation file) gives a pair:
the folder of the file in the config : <prefix>/<folder of the file in the design/reference folder>.
Values that already contain a placeholder such as {{REFERENCE_DIRECTORY}} are skipped.
The script exits with an error if an absolute path is not in the validation file, or if a file
has a different name in the validation file than in the config.

Output formats:
  pairs  one old:new pair per line (default)
  flags  one `--reference-path old:new` per line, for `hydra-genetics prepare-environment reference-path-update`
"""
import argparse
import posixpath
import sys

import yaml

PLACEHOLDER = "{{"
BEEGFS_PREFIX = "/beegfs-storage"
DEFAULT_PREFIX = "{{DESIGN_REF_FOLDER}}"


def walk_validation(node, keys=()):
    """Yield (key, entry) for every leaf of the validation file that has a url and a path."""
    if isinstance(node, dict):
        if "url" in node and "path" in node:
            yield "/".join(keys), node
        else:
            for key, value in node.items():
                yield from walk_validation(value, keys + (str(key),))


def walk_config(node, keys=(), in_list=False):
    """Yield (key, value, in_list) for every string in the config."""
    if isinstance(node, dict):
        for key, value in node.items():
            yield from walk_config(value, keys + (str(key),), in_list)
    elif isinstance(node, list):
        for index, value in enumerate(node):
            yield from walk_config(value, keys + (str(index),), True)
    elif isinstance(node, str):
        yield "/".join(keys), node, in_list


def url_to_path(url):
    """Turn file:/beegfs-storage/data/x into /data/x so it can be compared with config paths."""
    path = url[len("file:"):] if url.startswith("file:") else url
    if path.startswith(BEEGFS_PREFIX):
        path = path[len(BEEGFS_PREFIX):]
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--validation-file", required=True,
                        help="hydra-genetics references validation file (yaml)")
    parser.add_argument("-c", "--configfile", action="append", required=True,
                        help="config file to check, can be given more than once")
    parser.add_argument("-p", "--new-path-prefix", default=DEFAULT_PREFIX,
                        help=f"path the design/reference folder is replaced by (default: {DEFAULT_PREFIX})")
    parser.add_argument("-o", "--output", choices=["pairs", "flags"], default="pairs",
                        help="output output (default: pairs)")
    args = parser.parse_args()
    prefix = args.new_path_prefix.rstrip("/")

    with open(args.validation_file) as stream:
        validation = yaml.safe_load(stream)

    # config path -> design folder paths; a downloaded folder is its own entry, a file is in a folder
    downloads = {}
    folders = set()
    for _, entry in walk_validation(validation):
        url_path = url_to_path(entry["url"])
        if entry.get("type") == "folder":
            # a folder is downloaded from an archive inside it, so compare with the archive's folder
            url_path = posixpath.dirname(url_path)
            folders.add(url_path)
        downloads.setdefault(url_path, set()).add(entry["path"])

    pairs = {}
    errors = []
    for configfile in args.configfile:
        with open(configfile) as stream:
            config = yaml.safe_load(stream)
        for key, value, in_list in walk_config(config):
            if not value.startswith("/") or PLACEHOLDER in value:
                continue
            where = f"{configfile}\t{key}\t{value}"
            if value not in downloads:
                errors.append(f"{where}\tnot in the validation file")
                continue
            if len(downloads[value]) > 1:
                errors.append(f"{where}\tmatches several validation entries: {', '.join(sorted(downloads[value]))}")
                continue
            (design_path,) = downloads[value]
            if posixpath.basename(design_path) != posixpath.basename(value):
                errors.append(f"{where}\tis called {posixpath.basename(design_path)} in the validation file")
                continue
            if in_list:
                print(f"WARNING: {where}\tis inside a list, reference-path-update will not replace it", file=sys.stderr)
            if value in folders:
                old, new = value, f"{prefix}/{design_path}"
            else:
                design_dir = posixpath.dirname(design_path)
                old = posixpath.dirname(value)
                new = f"{prefix}/{design_dir}" if design_dir else prefix
            if pairs.setdefault(old, new) != new:
                errors.append(f"{where}\t{old} is mapped to both {pairs[old]} and {new}")

    if errors:
        print("Cannot create the path pairs:", file=sys.stderr)
        for error in errors:
            print(f"  {error}", file=sys.stderr)
        sys.exit(1)

    # longest first, since reference-path-update replaces substrings in the order given
    for old in sorted(pairs, key=lambda path: (-len(path), path)):
        pair = f"{old}:{pairs[old]}"
        print(f"--reference-path {pair}" if args.output == "flags" else pair)


if __name__ == "__main__":
    main()
