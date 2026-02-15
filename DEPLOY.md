# GitHub Pages Deployment Guide

This guide explains how to set up GitHub Pages for the scVIZ documentation.

## Prerequisites

- Repository pushed to GitHub at `https://github.com/jdwidawski/scVIZ`
- GitHub Actions enabled for the repository

## Setup Steps

### 1. Enable GitHub Pages

1. Go to your repository on GitHub: `https://github.com/jdwidawski/scVIZ`
2. Click on **Settings** (top navigation bar)
3. In the left sidebar, click on **Pages** (under "Code and automation")
4. Under "Build and deployment":
   - **Source**: Select "GitHub Actions"
   - This allows the `.github/workflows/docs.yml` workflow to deploy the docs

### 2. Verify Workflow Permissions

1. Still in **Settings**, go to **Actions** > **General** (left sidebar)
2. Scroll to "Workflow permissions"
3. Ensure "Read and write permissions" is selected
4. Check "Allow GitHub Actions to create and approve pull requests"
5. Click **Save**

### 3. Deploy Documentation

The documentation will automatically deploy when you push to the `master` branch:

```bash
cd /home/jakub_widawski/projects/scrnaseq-workflow
git push origin master
```

You can also manually trigger the workflow:
1. Go to **Actions** tab in your repository
2. Click on "Build and Deploy Documentation" workflow
3. Click "Run workflow" button
4. Select the `master` branch and click "Run workflow"

### 4. Access Your Documentation

After the workflow completes (usually 2-5 minutes):
- Documentation URL: `https://jdwidawski.github.io/scVIZ/`

### 5. Monitor Deployment

1. Go to the **Actions** tab in your repository
2. You should see the "Build and Deploy Documentation" workflow running
3. Click on the workflow run to see the build logs
4. Once complete, check the deployment status under **Settings** > **Pages**

## Workflow Details

The documentation deployment workflow (`.github/workflows/docs.yml`):

- **Triggers**: Automatically on push to `master` or `main` branch, or manual trigger
- **Process**:
  1. Checks out the code
  2. Sets up Python 3.11 and UV package manager
  3. Installs dependencies including Sphinx, RTD theme, and MyST parser
  4. Builds HTML documentation from `docs/source/` to `docs/build/html/`
  5. Creates `.nojekyll` file for GitHub Pages
  6. Uploads the HTML as a GitHub Pages artifact
  7. Deploys to GitHub Pages

## Troubleshooting

### Documentation not updating
- Check the workflow logs in the Actions tab
- Ensure the workflow completed successfully
- GitHub Pages may take a few minutes to update after deployment

### 404 Error
- Verify GitHub Pages is configured to use "GitHub Actions" as source
- Check that the workflow has "pages: write" permission
- Ensure `.nojekyll` file is present in the build output

### Broken links or missing pages
- Run the build locally: `cd docs && sphinx-build -b html source build/html`
- Check for warnings or errors in the build output
- Verify all `.rst` files are included in the table of contents

## Local Testing

To test the documentation locally before pushing:

```bash
# Activate virtual environment
source .venv/bin/activate

# Build documentation
cd docs
sphinx-build -b html source build/html

# Serve locally
cd build/html
python3 -m http.server 8080
```

Then open `http://localhost:8080` in your browser.

## Updating Documentation

To update the documentation:

1. Edit files in `docs/source/` (`.rst` or `.md` files)
2. Build locally to test: `sphinx-build -b html source build/html`
3. Commit changes: `git commit -m "docs: description of changes"`
4. Push to GitHub: `git push origin master`
5. GitHub Actions will automatically rebuild and deploy

## Badge Status

The README includes a documentation badge that shows build status:

```markdown
[![Documentation](https://github.com/jdwidawski/scVIZ/actions/workflows/docs.yml/badge.svg)](https://jdwidawski.github.io/scVIZ/)
```

This badge will be:
- 🟢 Green: Documentation built and deployed successfully
- 🔴 Red: Build failed (check workflow logs)
- 🟡 Yellow: Build in progress
