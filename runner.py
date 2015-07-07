import logging


logger=logging.getLogger(__file__)


class RunnerFSM(object):
    def __init__(self, dynamics, observer):
        self.dynamics=dynamics
        self.observer=observer

    def init(self):
        self.dynamics.init()

    def run(self):
        running=True
        while running:
            transition, when=self.dynamics.next()
            if transition is not None:
                self.dynamics.fire(transition, when)
                running=self.observer(transition, when)
            else:
                running=False
