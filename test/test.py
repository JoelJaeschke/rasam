import os
import io
import sys
import glob
import importlib
import importlib.util
import contextlib
import argparse

LEN_HEADER = 80
TEST_PATHS = os.path.join(os.path.dirname(__file__), "test_*.py")
TEST_FILES = glob.glob(TEST_PATHS)

class TestState:
    def __init__(self, **kwargs):
        self.state = kwargs

    def set(self, k, v):
        self.state[k] = v

    def get(self, k):
        if k in self.state:
            return self.state[k]
        else:
            return None

STATE = TestState()

class Colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def _module_from_path(path):
    file_name = path.split("/")[-1]
    module_name = file_name[:-3]

    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)

    return module_name, module

def _extract_tests_from_module(module):
    test_funcs = {
        key: value for key, value in module.__dict__.items() if key[:2] != "__" and key[:4] == "test"
    }

    return test_funcs

def _print_test_summary(summary):
    contains_failed_tests = False
    
    for i, (module, results) in enumerate(summary.items()):
        successful_funcs = {func: "success" for func in results["success"]}
        failed_funcs = {func: "fail" for func in results["failed"]}
        method_with_result = successful_funcs | failed_funcs
        
        module_successful = len(failed_funcs) == 0
        num_failed = len(failed_funcs)
        num_total = len(method_with_result)
        if module_successful:
            success_text = f"{Colors.OKGREEN}Success{Colors.ENDC}"
        else:
            success_text = f"{Colors.FAIL}Failure ({num_failed} of {num_total} failed){Colors.ENDC}"
            contains_failed_tests = True

        print(f"Results for module '{Colors.BOLD}{module}{Colors.ENDC}': {success_text}")
        for method, status in method_with_result.items():
            text = f"{Colors.OKGREEN}✔{Colors.ENDC}" if status == "success" else f"{Colors.FAIL}✗{Colors.ENDC}"
            print(f"\t{method}: {text}")

        if i < len(summary) - 1:
            print()

    if not contains_failed_tests:
        print()
        print(f"{Colors.OKGREEN}" + "="*LEN_HEADER)
        print(" "*27 + "All tests ran successfully!")
        print("="*LEN_HEADER + f"{Colors.ENDC}", end=None)
    else:
        print()
        print(f"{Colors.FAIL}" + "="*LEN_HEADER)
        print(" "*32 + "Some tests failed!")
        print("="*LEN_HEADER + f"{Colors.ENDC}", end=None)


def _print_test_io(name, stdout, stderr):
    stdout_content = stdout.rstrip()
    stderr_content = stderr.rstrip()

    if len(stdout_content) == 0 and len(stderr_content) == 0:
        return

    len_remaining = LEN_HEADER - len(name) - 2
    len_right = len_remaining // 2
    len_left = len_remaining - len_right
    header = "="*len_left + f" {name} " + "="*len_right
    
    print(header)
    
    print("Stdout:")
    if len(stdout_content) > 0:
        print(stdout_content)
    
    print("Stderr:")
    if len(stderr_content) > 0:
        print(stderr_content)
    
    print("="*LEN_HEADER)
    print()

def run_all_tests(paths, state, show_io = False, show_io_on_fail = False):
    test_results = {}

    encountered_failed_module = False    
    for path in paths:
        mod_name, module = _module_from_path(path)
        methods_to_run = _extract_tests_from_module(module)

        mod_result = {
            "success": [],
            "failed": []
        }

        for test, func in methods_to_run.items():
            failed = False
            with contextlib.redirect_stdout(io.StringIO()) as f_stdout, contextlib.redirect_stderr(io.StringIO()) as f_stderr:
                try:
                    func(state)
                    mod_result["success"].append(test)
                except:
                    mod_result["failed"].append(test)
                    failed = True
                    encountered_failed_module = True

            if (failed and show_io_on_fail) or show_io:
                _print_test_io(test, f_stdout.getvalue(), f_stderr.getvalue())

        test_results[mod_name] = mod_result

    _print_test_summary(test_results)

    if encountered_failed_module:
        return 1
    else:
        return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog='Rasam test driver',
                    description='Run all tests for rasam without requiring any imports')
    parser.add_argument('-v', dest="verbose", action="store_true")
    parser.add_argument('-vv', dest="more_verbose", action="store_true")
    parser.add_argument('-p', '--path', dest="tool_path", help="Path to rasam tool. If installed, just use 'rasam', if built from scratch, point to the binary.")
    args = parser.parse_args()

    if args.tool_path:
        STATE.set("UTIL_PATH", args.tool_path)
    else:
        STATE.set("UTIL_PATH", "/host/build/src/rasam")

    sys.exit(run_all_tests(TEST_FILES, STATE, show_io=args.more_verbose, show_io_on_fail=args.verbose))