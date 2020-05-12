class BaseService(object):
    def __init__(self):
        pass

    def process_args(self, args):
        raise NotImplementedError()

    def run(self):
        raise NotImplementedError()
