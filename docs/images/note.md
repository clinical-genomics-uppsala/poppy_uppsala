# Best Practice: Place Images Inside the docs/ Directory

MkDocs only collects static files (such as images) that are inside the docs/ directory.
Alternatively, you can use a static link to the raw.githubusercontent URL of the image, but this is not recommended for 
several reasons:

1. Performance: Serving images from your repository can be slower than serving them from a CDN.
2. Availability: If the repository is private or the image is deleted, it will not be accessible.
3. Consistency: Using local images ensures that they are always available and consistent across different environments.
4. Documentation Clarity: Keeping images within the docs/ directory makes it clear that they are part of the documentation and not just random files in the repository.

Example of a static link to an image:
`![rulegraph](https://raw.githubusercontent.com/clinical-genomics-uppsala/poppy_uppsala/patch-readthedocs/images
/rulegraph.png)`

Example of a local image in the docs/ directory:
`![rulegraph](images/rulegraph.png)`

If your images are outside root/docs/ (for example, at the repository root or in a subdirectory like root/images/),
they will not be included in the site build and
therefore will not display on Read the Docs or in your deployed MkDocs site.

This is highlighted in discussions like @squidfunk/mkdocs-material/discussions/3424.

Recommended folder structure:
```
your-repo/
├── docs/
│   ├── index.md
│   ├── usage.md
│   └── images/
│       └── your-image.png
├── mkdocs.yaml
└── ...
```
