# Contributing

## Reporting bugs 🐛

Before reporting a bug, try to search for a similar problem in the *Issues* section on GitHub. Clear the search bar to include closed issues in the search and type your phrase.

If you have not found a relevant issue, you can report the bug in the *Issues* section. When you click on *New issue*, various templates will be displayed — please pick *Bug Report*. Carefully fill out the template and submit the issue.

## Requesting features 💡

If you have an idea for improvement, you can submit your proposal in the *Issues* section on GitHub. When you click on *New issue*, various templates will be displayed — please pick *Feature Request*. Carefully fill out the template and submit the issue.

## Pull requests ⤴️

**Working on your first Pull Request?** You can learn how from this *free* series [How to Contribute to an Open Source Project on GitHub](https://kcd.im/pull-request)

### Local development environment setup

In order to develop/test the library one should have the following dependencies
installed:
* _g++/clang++_
* _cpplint_
* _doxygen_
* _make_
* _mpfr_
* _pre-commit_
* _valgrind_

For Ubuntu we recommend `apt` package manager, for macOS - `brew`.

Following cloning the repository the first action one should take is to
install pre-commit hooks for `git` with the following command:
```bash
pre-commit install
```
Dev commands have been encapsulated in the _Makefile_, see: `make help` for more info.

An alternative approach to set up the environment involves building
an isolated environment, just for this project:

1. Install Conda (package management system and environment management). You can find the installation instructions in
   [Anaconda page](https://www.anaconda.com/) or in [Miniconda page](https://docs.conda.io/en/latest/miniconda.html).
2. Setup a Conda environment to install all needed dependencies:
   ```bash
   conda env create
   ```
3. In order to activate newly created Conda environment type:

   ```bash
   conda activate hypercomplex-dev
   ```

   If you have no experience with Conda environments, read about them
   [under this address](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html).

### local Docker container

We have additionally prepared a development/execution Docker image which one may use in order to
work on our project in a fully encapsulated environment (that is: a container).  
Assuming the Docker Engine is running locally please build the image with:
```
docker build -t hypercomplex:latest -f Dockerfile-dev .
```
Once the process is finished start & enter the container with:
```
docker run --name hypercomplex -e HOSTUID=`id -u $USER` -p 8888:8888 -v {HOST_PATH_TO_HYPERCOMPLEX}:/hypercomplex -it fuzzyreg:latest
```
Watch out! Due to Docker's specifics they need to be executed as `root` user;
[alternatively, see here](https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user).

Recall that all data generated inside the container (with the exception of the mounted volume) are **not** persistent.  
If you'd like your data don't perish into oblivion after you stop the container
check out [Docker documentation on storage mechanisms](https://docs.docker.com/storage/).

### Ephemeral development environment

(╯°□°)╯︵ ┻━┻  
If that above is you, feeling like flipping a table just by the sheer looks of all the instructions above, do not despair! You still can contribute to our codebase from the cloud! We set up an ephemeral [Gitpod](https://www.gitpod.io) environment for all the developers who prefer coding from a remote server.

Just click on this cool button below:  
[![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/AngryMaciek/hypercomplex)

See more information on how to start up a _Gitpod_ environment dedicated to a specific remote branch [here](https://www.gitpod.io/docs/introduction/learn-gitpod/context-url#branch-and-commit-contexts), specific issue [here](https://www.gitpod.io/docs/introduction/learn-gitpod/context-url#issue-context) and a specific pull request [here](https://www.gitpod.io/docs/introduction/learn-gitpod/context-url#pullmerge-request-context).

However, please remember that such luxury is limitted:

> Gitpod offers a free plan for new users which includes 50 hours of standard workspace usage.
> If you need more hours, you can upgrade to one of the paid plans in your personal settings.

### Branch naming convention

You can read about git branching [under this address](https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging). The branch names should follow the convention specified below:

```
<type>/<issue id>/<short description>
```

where

- **type** is one of the [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/)' types, i.e. `feat`, `ci`, `docs`, `fix`, `refactor`, `test`, `chore`, etc.
- **issue id** is the id of the issue which will be fixed by this PR.
- **short description** is ~25 characters long at max and is hyphenated.

Examples:

```
chore/6/gitignore
ci/20/testing-and-code-coverage
ci/21/flake8-action
ci/46/docstring-linting
docs/22/contributing-instructions
feat/3/license
feat/4/github-templates
feat/40/version-flag
feat/7/code-of-conduct
refactor/8/initial-project-structure
```

### Code style & testing

To make it easier for everyone to maintain, read and contribute to the code,
as well as to ensure that the code base is robust and of high quality, we
would kindly ask you to stick to the guidelines for code style and
testing. GitHub Actions mechanism of Continous Integration tests the code
style with _cpplint_ and performs unit tests with the _Catch2_ framework
(measuring code coverage alongside).
Also, please try to conform to the used code, docstring and commenting style within
a project to maintain consistency.

### Merging the pull request

A pull request can only be merged after it completed all of the checks performed by CI. After a manual review by one of the maintainers, it can be merged into `master`. The merge commit should follow the rules outlined in [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/).
