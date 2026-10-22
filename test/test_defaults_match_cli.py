"""A command function called with only its main input behaves like the command given only its
input (and an output): every keyword named like one of the command's options defaults to that
option's default, and the options the command requires are required by the function too.
"""
import importlib
import inspect

import pytest

from domainator.cli.registry import TOOLS

# Options the command requires that the function can do without.
OPTIONAL_IN_THE_FUNCTION = {
    ("hmmer_build", "name"): "None keeps the name the MSA file gives the alignment",
    ("build_ssn", "xgmml_out"): "the function can write only the cluster table",
}


@pytest.mark.parametrize("name", list(TOOLS))
def test_keyword_defaults_match_the_options(name):
    module = importlib.import_module(f"domainator.{TOOLS[name].module}")
    parser = module.build_parser()
    actions = {action.dest: action for action in parser._actions}
    declared_outputs = set(getattr(parser, "required_outputs", None) or ())
    for parameter in inspect.signature(getattr(module, name)).parameters.values():
        action = actions.get(parameter.name)
        if parameter.kind is not inspect.Parameter.KEYWORD_ONLY or action is None:
            continue
        has_default = parameter.default is not inspect.Parameter.empty
        if action.required:
            assert not has_default or (name, parameter.name) in OPTIONAL_IN_THE_FUNCTION, parameter.name
        elif not has_default:
            # the one output a command must write can be required by its function
            assert parameter.name in declared_outputs, f"{parameter.name} has no default (the option's is {action.default!r})"
        else:
            assert parameter.default == action.default, \
                f"{parameter.name}={parameter.default!r}, but the option's default is {action.default!r}"
