import difflib
import copy
import itertools
import yaml
import os
import gzip
import pybedtools
from .tfuncs import setup_module, teardown_module

yamltestdesc = ["test_cases.yaml"]

if pybedtools.settings._v_2_27_plus:
    yamltestdesc.append("test_merge227.yaml")
    yamltestdesc.append("test_shuffle227.yaml")

elif pybedtools.settings._v_2_15_plus and not pybedtools.settings._v_2_27_plus:
    yamltestdesc.append("test_merge215.yaml")
    yamltestdesc.append("test_shuffle215.yaml")

this_dir = os.path.dirname(__file__)
yamltestdesc = [os.path.join(this_dir, i) for i in yamltestdesc]


def gz(x):
    """
    Gzips a file to a tempfile, and returns a new BedTool using the gzipped
    version.
    """
    gzfn = pybedtools.BedTool._tmp()
    with gzip.open(gzfn, "wb") as out_:
        with open(x.fn, "rb") as in_:
            out_.writelines(in_)
    return pybedtools.BedTool(gzfn)


def fix(x):
    """
    Replaces spaces with tabs, removes spurious newlines, and lstrip()s each
    line. Makes it really easy to create BED files on the fly for testing and
    checking.
    """
    s = ""
    for i in x.splitlines():
        i = i.strip("\n\r")
        if len(i) == 0:
            continue

        # If the expected output contains tabs, then use those to split,
        # otherwise space.  This allows you to have expected output with blank
        # fields (e.g., "\t\t")
        if "\t" in i:
            i = i.split("\t")
        else:
            i = i.split()

        i = "\t".join(i) + "\n"
        s += i
    return s


# List of methods that *only* take BAM as input
bam_methods = ("bam_to_bed",)

# List of supported BedTool construction from BAM files.  Currently only
# file-based.
supported_bam = ("filename",)

converters = {
    "filename": lambda x: pybedtools.BedTool(x.fn),
    "generator": lambda x: pybedtools.BedTool(i for i in x),
    "stream": lambda x: pybedtools.BedTool(open(x.fn)),
    "gzip": gz,
}


def run(d):
    method = d["method"]
    bedtool = d["bedtool"]
    convert = d["convert"]
    kwargs = d["kw"].copy()
    expected = d["test_case"]["expected"]

    bedtool_converter = convert.pop("bedtool")
    bedtool = converters[bedtool_converter](pybedtools.example_bedtool(bedtool))

    for k, converter_name in convert.items():
        kwargs[k] = converters[converter_name](pybedtools.example_bedtool(kwargs[k]))
    result = getattr(bedtool, method)(**kwargs)
    res = str(result)
    expected = fix(expected)
    try:
        assert res == expected

    except AssertionError:
        print(result.fn)
        print("Method call:")
        args = []
        for key, val in list(kwargs.items()):
            args.append(("%s=%s" % (key, val)).strip())

        args = ", ".join(args)
        print("BedTool.%(method)s(%(args)s)" % locals())
        print("Got:")
        print(res)
        print("Expected:")
        print(expected)
        print("Diff:")
        for i in difflib.unified_diff(res.splitlines(1), expected.splitlines(1)):
            print(i, end=" ")

        # Make tabs and newlines visible
        spec_res = res.replace("\t", "\\t").replace("\n", "\\n\n")
        spec_expected = expected.replace("\t", "\\t").replace("\n", "\\n\n")

        print("Showing special characters:")
        print("Got:")
        print(spec_res)
        print("Expected:")
        print(spec_expected)
        print("Diff:")
        for i in difflib.unified_diff(
            spec_res.splitlines(1), spec_expected.splitlines(1)
        ):
            print(i, end=" ")
        raise


def pytest_generate_tests(metafunc):
    tests = []
    labels = []
    for config_fn in yamltestdesc:
        if hasattr(yaml, "FullLoader"):
            test_cases = yaml.load(open(config_fn).read(), Loader=yaml.FullLoader)
        else:
            test_cases = yaml.load(open(config_fn).read())
        for test_case in test_cases:
            kw = test_case["kwargs"]

            kwc = copy.copy(kw)

            method = test_case["method"]

            a_isbam = False
            b_isbam = False
            i_isbam = False

            # Figure out if this is a test involving a method that operates on
            # two files (ab), one file (i) or the "bed" tests (bed)
            flavor = None
            if (("a" in kw) and ("b" in kw)) or ("abam" in kw):
                flavor = "ab"
            if ("i" in kw) or ("ibam" in kw):
                flavor = "i"

            # If bams were specified in the test block we need to keep track of
            # that. Then we set the 'a' or 'i' kws, which control the nature of
            # the constructed BedTool, to that filename.
            if "abam" in kw:
                kw["abam"] = pybedtools.example_filename(kw["abam"])
                kw["a"] = kw["abam"]
                a_isbam = True

            if "ibam" in kw:
                kw["ibam"] = pybedtools.example_filename(kw["ibam"])
                kw["i"] = kw["ibam"]
                i_isbam = True

            kinds = ["filename", "generator", "stream", "gzip"]

            if "files" in kw:
                kw["files"] = [pybedtools.example_filename(i) for i in kw["files"]]

            if "bams" in kw:
                kw["bams"] = [pybedtools.example_filename(i) for i in kw["bams"]]

            if "fi" in kw:
                kw["fi"] = pybedtools.example_filename(kw["fi"])

            if flavor == "i":
                orig_i = pybedtools.example_bedtool(kw["i"])
                if orig_i._isbam:
                    i_isbam = True
                bedtool = kw.pop("i")
                for kind in kinds:
                    if i_isbam and (kind not in supported_bam):
                        continue
                    label = "{method}: {kwc} {kind}".format(**locals())
                    labels.append(label)
                    tests.append(
                        dict(
                            method=method,
                            bedtool=bedtool,
                            test_case=test_case,
                            kw=kw,
                            convert={"bedtool": kind},
                        )
                    )

            if flavor == "ab":
                orig_a = pybedtools.example_bedtool(kw["a"])
                orig_b = pybedtools.example_bedtool(kw["b"])
                if orig_a._isbam:
                    a_isbam = True
                if orig_b._isbam:
                    b_isbam = True

                bedtool = kw.pop("a")

                for kind_a, kind_b in itertools.permutations(kinds, 2):
                    label = "{method}: {kwc} a={kind_a} b={kind_b}".format(**locals())
                    if a_isbam and (kind_a not in supported_bam):
                        continue
                    if b_isbam and (kind_b not in supported_bam):
                        continue

                    labels.append(label)
                    tests.append(
                        dict(
                            method=method,
                            bedtool=bedtool,
                            test_case=test_case,
                            kw=kw,
                            convert={"bedtool": kind_a, "b": kind_b},
                        )
                    )

    metafunc.parametrize("tests", tests, ids=labels)


def test_all(tests):
    run(tests)
