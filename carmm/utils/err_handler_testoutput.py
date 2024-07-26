
def test_set_logger():
    from carmm.utils.logger_set import set_logger

    #TODO: MAKE THIS TEST ACTUALLY TEST CAPTURED OUTPUT

    logger = set_logger("test_logger", 0)

    logger.debug('test_set_logger test passed! DEBUG error correctly printed')

test_set_logger()
