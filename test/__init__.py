import logging
import tempfile
import os

TEST_LOG = os.path.join(tempfile.gettempdir(), 'sagar_unittest.log')
try:
    os.remove(TEST_LOG)
except OSError:
    pass
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)s()] %(message)s"
logging.basicConfig(filename=TEST_LOG, level=logging.INFO, format=FORMAT)
