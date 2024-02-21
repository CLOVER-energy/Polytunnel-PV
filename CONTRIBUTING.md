# Contributing to Polytunnel-PV

:+1::tada: First off, thanks for taking the time to contribute! :tada::+1:

The following is a set of guidelines for contributing to Polytunnel-PV, which are hosted in the [CLOVER-energy Organization](https://github.com/CLOVER-energy) on GitHub. These are mostly guidelines, not rules. Use your best judgment, and feel free to propose changes to this document in a pull request.

#### Table Of Contents

[Code of Conduct](#code-of-conduct)

* [What should I know before I get started?](#what-should-i-know-before-i-get-started)
  * [What is Polytunnel-PV for?](#what-is-polytunnel-pv-for)
  * [What is Polytunnel-PV not for?](#what-is-polytunnel-pv-not-for)

* [How to contribute to Polytunnel-PV](#how-to-contribute-to-polytunnel-pv)
  * [Reporting bugs](#reporting-bugs)
  * [Merging patches](#merging-patches)
    * [Cosmetic patches](#cosmetic-patches)
  * [Questions](#questions)
  * [Contributing to the documentation](#contributing-to-the-documentation)

* [Styleguides](#styleguides)
  * [Git commit messages](#git-commit-messages)
  * [Python styleguide](#python-styleguide)
  * [YAML styleguide](#yaml-styleguide)
  * [Changing styleguides](#changing-styleguides)

* [Additional Notes](#additional-notes)
  * [Issue and pull request labels](#issue-and-pull-request-labels)


### Code of Conduct

This project and everyone participating in it is governed by the [Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code. Please report unacceptable behavior to the CLOVER development team.

### What should I know before I get started?

Polytunnel-PV was developed at Imperial College London as a means of simulating the performance of curved photovoltaic 
(PV) panels, primarily for applications in agricultural-photovoltaic (agri-PV) installations. Under continuous development since 2023, Polytunnel-PV is continually changing and improving to enable a greater number of integrations and applications to be considered.

#### What is Polytunnel-PV for?

Polytunnel-PV is for simulating the performance of curved PV modules. It is for assessing the impact of bypass diodes (diodes which enable current to flow past cells which are partially or fully shaded and so are producing less power than those under full illumination) and for determining the optimum configuration of these systems.

#### What is Polytunnel-PV not for?

Polytunnel-PV looks at the performance of a single, isolated solar module. It is not for the simulation of wider energy systems (though it can be integrated with the [CLOVER](https://github.com/CLOVER-energy/CLOVER) energy-system model) and so any performance calculations are those of the system operating under optimal maximum-power-point (MPP) conditions. The results can be used to help aid the _design_ of these modules, and give a rough indication of their performance, but should not be used for accurate sizing of systems.

## How to contribute to Polytunnel-PV

### Reporting bugs

**Did you find a bug?** Bugs make it into our code from time to time. If you spot a bug, report it as follows:

* **Ensure the bug was not already reported** by searching on GitHub under [Issues](https://github.com/CLOVER-energy/Polytunnel-PV/issues).

* If you're unable to find an open issue addressing the problem, [open a new one](https://github.com/CLOVER-energy/Polytunnel-PV/issues/new/choose). Be sure to include a **title and clear description**, as much relevant information as possible, and a **code sample** or an **executable test case** demonstrating the expected behavior that is not occurring.

  * If the issue is a **bug**, use the [Bug report](https://github.com/CLOVER-energy/Polytunnel-PV/issues/new?assignees=&labels=bug&template=bug_report.md&title=) template,

  * If the issue is a **feature** request for something new that you would like to see introduced into Polytunnel-PV, use the [Feature request](https://github.com/CLOVER-energy/Polytunnel-PV/issues/new?assignees=&labels=enhancement&template=feature_request.md&title=) template.

### Merging patches

**Did you write a patch that fixes a bug?** If you have coded a solution for a bug that you have found or for an open issue, open a pull request for it as follows:

* Open a new GitHub pull request with the patch.

* Ensure the PR description clearly describes the problem and solution:

  * Include the relevant issue number if applicable,

  * Follow the template information presented, filling in all the fields requested which are relevant to your patch.

* Ensure that you include at least one administrator reviewer for your pull request. Without an appropriate review, you will be unable to merge your pull request.

#### Cosmetic patches

**Did you fix whitespace, format code, or make a purely cosmetic patch?** Changes that are cosmetic in nature and do not add anything substantial to the stability, functionality, or testability of Polytunnel-PV will generally not be accepted. Contact the developers directly, or save your cosmetic changes until you are able to merge them as part of a feature or bug issue.

### Questions

**Do you have questions about the source code?** Ask any question about how to use Polytunnel-PV on the [Discussions](https://github.com/CLOVER-energy/Polytunnel-PV/discussions) page.

### Contributing to the documentation

**Do you want to contribute to the Polytunnel-PV documentation?** Polytunnel-PV is an ever-evolving piece of software. If you want to contribute to the documentation, get in touch with the development team. Documentation updates are usually produced for major releases.

## Styleguides

### Git commit messages

* Git Commit Messages
* Use the present tense ("Add feature" not "Added feature")
* Use the imperative mood ("Move cursor to..." not "Moves cursor to...")
* Limit the first line to 72 characters or less
* Reference issues and pull requests liberally after the first line
* When only changing documentation, include [ci skip] in the commit title
* Consider starting the commit message with an applicable emoji:
  * 🎉 `:tada` when adding an initial commit or solving a difficult problem
  * 🎨 `:art:` when improving the format/structure of the code
  * 🐎 `:racehorse:` when improving performance
  * 📝 `:memo:` when writing docs
  * 🐧 `:penguin:` when fixing something on Linux
  * ✨ `:sparkles:` when adding a new feature
  * 🚧 `:construction:` when part-way through working on code
  * 🍎 `:apple:` when fixing something on macOS
  * 🏁 `:checkered_flag:` when fixing something on Windows
  * ⏪ `:rewind:` when backing out a commit or changes
  * 🐛 `:bug:` when fixing a bug
  * 🔥 `:fire:` when removing code or files
  * 💚 `:green_heart:` when fixing the CI build
  * 👕 `:shirt:` when removing linter warnings
  * ✅ `:white_check_mark:` when adding tests
  * 🔒 `:lock:` when dealing with security
  * ⬆️ `:arrow_up:` when upgrading dependencies
  * ⬇️ `:arrow_down:` when downgrading dependencies
  * 🚀 `:rocket:` when deploying code

### Python styleguide

All `Python` code must conform to [mypy](https://github.com/python/mypy) and [pylint](https://github.com/PyCQA/pylint) coding standards and must be formatted with the [black](https://github.com/psf/black) formatter:
* A `mypy.ini` file within the root of the repository sets out the requirements that must be met by any code. Use `python -m mypy src/` to ensure that your code complies with the guidelines.
* A `.pylintrc` file within the root of the repository sets out the linting requirements. Use `python -m pylint src/` to ensure that your code complies with the guidelines.
* All code must be formatted with `black`.

These tests must pass for any pull request to be successfully approved and merged. You can run these tests from the root of the repository with `./bin/test-polytunnel-pv.sh`.

### YAML styleguide

All `.yaml` files which are modified are linted with [yamllint](https://github.com/adrienverge/yamllint). You can use `yamllint -c .yamllint-config.yaml` to run `yamllint` over any `.yaml` files that you have modified.

### Changing styleguides

If you have any changes for the styleguides, make these **very clear** within your pull request message. It is unusual to have to change the styleguides to make your code pass the required tests.

## Additional notes

Thanks! :heart: :heart: :heart:

CLOVER-energy Development Team
