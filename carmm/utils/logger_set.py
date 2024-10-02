def set_logger(logger_name, warning_level = 1):

    """
    A general purpose error handler for CARMM - sets the Logger object
    with a sensible set of defaults
    
    Args:
        logger_name: str
            Name of logger passed to get_logger
        warning_level: int
            Sets minimum-severity log message printed from lowest (0, DEBUG) to highest (4, CRITICAL)

    Returns:
        carmm_logger: logger
            Logger object to be used by other modules for error handling

    """ 

    import logging
    from sys import stdout

    # Create dictionary for easy access debugs
    warning_lvls = {0: logging.DEBUG, 1: logging.INFO, 2: logging.WARNING,
                    3: logging.ERROR, 4: logging.CRITICAL}

    # create logger
    logger = logging.getLogger(logger_name)
    logger.setLevel(warning_lvls[warning_level])

    # create console handler and set level to debug
    ch = logging.StreamHandler(stdout)
    ch.setLevel(warning_lvls[warning_level])

    # create formatter
    # FULL FORMATTING
    #formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # REDUCED FORMATTING
    formatter = logging.Formatter('%(levelname)s - %(message)s')

    # add formatter to ch
    ch.setFormatter(formatter)

    # add ch to logger
    logger.addHandler(ch)

    return logger

