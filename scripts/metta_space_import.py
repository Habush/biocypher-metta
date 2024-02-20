from hyperon import *
import typer
from typing_extensions import Annotated
import pathlib
import os
import time
import datetime
import resource
import logging

app = typer.Typer()

class Timer(object):
    def __init__(self, name=None, logger_name=None):
        self.name = name
        self.logger = logging.getLogger(logger_name)

    def __enter__(self):
        self.tstart = time.time()

    def __exit__(self, type, value, traceback):
        message = 'Elapsed time: %s' % (datetime.timedelta(seconds=time.time() - self.tstart))
        if self.name:
            message = '[%s] ' % self.name + message
        self.logger.debug(message)


def memory_usage(point=""):
    usage=resource.getrusage(resource.RUSAGE_SELF)
    return '''%s: usertime=%.3f systime=%.3f mem=%.3f mb
           '''%(point,usage[0],usage[1],
                usage[2]/1024.0 )


def setup_logger(logger_name, log_file, level=logging.INFO):
    logger = logging.getLogger(logger_name)
    formatter = logging.Formatter('%(asctime)s - %(name)s : %(message)s')
    fileHandler = logging.FileHandler(log_file, mode='w')
    fileHandler.setFormatter(formatter)
    streamHandler = logging.StreamHandler()
    streamHandler.setFormatter(formatter)

    logger.setLevel(level)
    logger.addHandler(fileHandler)
    logger.addHandler(streamHandler)

    return logger

@app.command()
def load_metta_space(input_dir: Annotated[pathlib.Path,
                        typer.Option(exists=True, file_okay=False, dir_okay=True)],
                     type_def_path: Annotated[pathlib.Path,
                        typer.Option(exists=True, file_okay=True, dir_okay=False)],
                     log = None):

    if log and os.path.exists(log):
        os.remove(log)

    logger = setup_logger('metta_space_import', log, logging.DEBUG)

    with Timer("Loading MeTTa space...", logger_name="metta_space_import"):
        metta = MeTTa(env_builder=Environment.test_env())
        logger.info(f"Loading type definitions ...")
        metta.import_file(str(type_def_path.resolve()))
        logger.debug(memory_usage("After loading type definitions"))
        for path in input_dir.rglob("*.metta"):
            full_path = str(path.resolve())
            logger.info(f"Loading {full_path} ...")
            metta.import_file(full_path)
            logger.debug(memory_usage(f"After loading {full_path}"))

        # get properties of (gene ENSG00000290825)
        prog1 = '''
            !(match &self ($x (gene ENSG00000290825) $y) ($x (gene ENSG00000290825) $y))
        '''
        with Timer(f"Executing query : {prog1}", logger_name="metta_space_import"):
            logger.info(metta.run(prog1))

        logger.debug(memory_usage("After executing query"))

        # find genes on chr16 b/n base numbers 53MB and 56MB
        prog2 = '''
            !(match &self (, (chr $g "chr16")
                            (start $g $start)
                            (end $g $end))
                    (if (and (> $start 53000000) (< $end 56000000)) $g ()))
        '''

        with Timer(f"Executing query : {prog2}", logger_name="metta_space_import"):
            logger.info(metta.run(prog2))

        logger.debug(memory_usage("After executing query"))


if __name__ == "__main__":
    app()