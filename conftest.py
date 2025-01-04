#Custom pass/fail message for pytest (because I'm extra)
import sys
try:
    import pytest
    def pytest_terminal_summary(terminalreporter, exitstatus, config):
        if exitstatus == 0:
            terminalreporter.write("\n🎉 ALL TESTS PASS 🎉\n") 
        else:
            terminalreporter.write("\n💥 SOME TESTS FAILED 💥\n")
            
except ImportError:
    print("pytest not found. Falling back to unittest.")
    
    def run_unittest():
        import unittest
        loader = unittest.TestLoader()
        suite = loader.discover('.')
        runner = unittest.TextTestRunner(verbosity=2)
        result = runner.run(suite)

        if result.wasSuccessful():
            print("\n🎉 ALL TESTS PASS 🎉")
        else:
            print("\n💥 SOME TESTS FAILED 💥")
        sys.exit(0 if result.wasSuccessful() else 1)

    if __name__ == "__main__":
        run_unittest()

