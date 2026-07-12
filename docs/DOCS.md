## Documentation
The documentation is built using Sphinx. Assuming you have Python 3 and [pipenv](https://pipenv.readthedocs.io/en/latest/install/#installing-pipenv)
installed, run the following:

```shell
cd docs
pipenv install --dev
```

#### Running the Development Server
You can run a local development server to preview changes using the following:

```shell
cd docs
pipenv run serve

```

#### Running Manually
Build with the warning and reference checks used in CI:

```shell
cd docs
pipenv run make strict
```

Check generated command and LLM artifacts and run the retrieval benchmark:

```shell
cd docs
pipenv run make check-generated
pipenv run make benchmark
```

Run one isolated command from every tutorial:

```shell
cd docs
pipenv run make smoke
```
