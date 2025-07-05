import argparse
from copy import deepcopy as dc
import json
import os
import re
import textwrap

from CABS.io import logger

_HELPW = 100
_wrapper = textwrap.TextWrapper(width=_HELPW, break_long_words=True, expand_tabs=False)


def _wrap(text):
    lines = []
    for line in text.split("\n"):
        _list = _wrapper.wrap(line)
        if not _list:
            _list = ["\n"]
        lines.extend(_list)
    return "\n".join(lines)


class CABSFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog, indent_increment=2, max_help_position=4, width=_HELPW):
        super(CABSFormatter, self).__init__(
            prog, indent_increment, max_help_position, width
        )

    def _split_lines(self, text, width):
        return super(CABSFormatter, self)._split_lines(text, width) + [" "]


class ConfigFileParser:
    OPTIONRE = re.compile(r"(?P<option>[^:=]*)" r"[:=]" r"(?P<value>.*)$")

    def __init__(self, filename):
        self.args = []
        with open(filename) as f:
            for line in f:
                if line == "" or line[0] in ";#\n":
                    continue
                match = self.OPTIONRE.match(line)
                if match:
                    option, value = match.groups()
                    self.args.append("--" + option.strip())
                    self.args.extend(value.split("#")[0].split(";")[0].split())
                else:
                    try:
                        test = ["--" + line.split("#")[0].split(";")[0].split()[0]]
                        self.args.extend(test)
                    except ValueError:
                        pass


def mk_usage(parser_dict, option_dict, indent=4):
    prog = parser_dict["prog"]
    usage = [f"usage: {prog} [OPTIONS]", ""]
    for group, _options in parser_dict["groups"]:
        usage.append(group)
        for name in _options:
            option = option_dict[name]
            flag = option.get("flag")
            metavar = option.get("metavar")
            line = " " * indent
            if flag:
                line += flag + ", "
            else:
                line += " " * 4
            line += "--" + name
            if metavar:
                if type(metavar) is tuple:
                    metavar = " ".join(metavar)
                line += " " + str(metavar)
            usage.append(line)
    usage.extend(["", f"For full help run: {prog} -h, --help"])
    return "\n".join(usage)


def mk_parser(parser_dict, group_dict, option_dict):
    parser_dict["description"] = _wrap(parser_dict["description"])
    parser_dict["epilog"] = _wrap(parser_dict["epilog"])

    _groups = parser_dict.pop("groups")
    defaults = parser_dict.pop("defaults", {})
    parser = argparse.ArgumentParser(
        formatter_class=CABSFormatter,
        add_help=False,
        usage="%s [OPTIONS]" % parser_dict["prog"],
        **parser_dict,
    )
    for group_name, _options in _groups:
        group = group_dict[group_name]
        group["description"] = _wrap(group["description"])
        group = parser.add_argument_group(title=group_name, **group)
        for opt_name in _options:
            option = option_dict[opt_name]
            name = ["--" + opt_name]
            flag = option.pop("flag", None)
            if flag:
                name.insert(0, flag)
            if opt_name in defaults:
                option["default"] = defaults[opt_name]
            option["help"] = _wrap(option["help"])
            group.add_argument(*name, **option)
    return parser


def restore_types(dict):
    from CABS.config_loader import get_type_dispatch

    TYPE_DISPATCH = get_type_dispatch()

    # Convert string type names to actual functions
    type_functions = {
        "int": int,
        "float": float,
        "str": str,
        "split_equals": lambda x: x.split("="),
    }

    for opt in dict.values():
        if "type" in opt and opt["type"] in type_functions:
            opt["type"] = type_functions[opt["type"]]
        if "type" in opt and isinstance(opt["type"], str):
            type_key = opt["type"]
            if type_key in TYPE_DISPATCH:
                opt["type"] = TYPE_DISPATCH[type_key]


def load_config(filename="config.json"):
    base_dir = os.path.dirname(__file__)
    config_path = os.path.join(base_dir, filename)

    with open(config_path) as f:
        config = json.load(f)

    restore_types(config["options"])

    return (
        config["dock_dict"],
        config["flex_dict"],
        config["options"],
        config["groups"],
    )


dock_dict, flex_dict, options, groups = load_config()

dock_usage = mk_usage(dock_dict, options)
flex_usage = mk_usage(flex_dict, options)
dock_parser = mk_parser(dock_dict, dc(groups), dc(options))
flex_parser = mk_parser(flex_dict, dc(groups), dc(options))


def if_append(option_name, value):
    """Handles appended arguments that come as both lists of lists and lists of arguments"""
    try:
        if options[option_name]["action"] == "append":
            try:
                nargs = options[option_name].get("nargs")
                if type(value) == list:
                    line = ""
                    for single_value in value:
                        line += (
                            "\n"
                            + option_name
                            + " : "
                            + " ".join([str(i) for i in single_value])
                        )
                    return line
                else:
                    logger.warning(
                        "OptParse",
                        "Issues while saving multiple argument option: %s"
                        % option_name,
                    )
                    raise KeyError

            except KeyError:
                if type(value) == list:
                    line = ""
                    for single_value in value:
                        line += "\n" + option_name + " : " + str(single_value)
                    return line
                else:
                    logger.warning(
                        "OptParse",
                        "Issues while saving appended argument option: %s"
                        % option_name,
                    )
                    raise KeyError
            except Exception:
                logger.warning(
                    "OptParse", "Issues while saving %s option" % option_name
                )
                raise KeyError
        else:
            raise KeyError
    except KeyError:
        raise


def if_store_true(option_name, value):
    """Catches flags that are not True"""
    try:
        if options[option_name]["action"] == "store_true":
            if value:
                return "\n" + option_name
            else:
                return " "
    except KeyError:
        raise


def if_nargs(option_name, value):
    """Handles options that come as lists"""
    try:
        nargs = options[option_name]["nargs"]
        return (
            "\n" + option_name + " : " + " ".join([str(i).strip("#") for i in value])
        )  # "#" for contact maps colors
    except KeyError:
        raise


def if_wd(option_name, value):
    if option_name == "work-dir":
        return "\n" + option_name + " : " + os.path.abspath(value)
    else:
        raise KeyError


special_cases = [if_append, if_store_true, if_nargs, if_wd]


def option_formatter(option, value):
    """
    Provides a string with properly formatted option to save in config file
    If value is False the option is ignored, this should be in line with options design
    """

    if value is None or not value:
        return " "
    for func in special_cases:
        try:
            return func(option, value)
        except KeyError:
            pass
    return "\n" + option + " : " + str(value)
