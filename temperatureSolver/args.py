import argparse

class Args():
    def __init__(self):

        parser = argparse.ArgumentParser()
        parser.add_argument('--height', type=float, default=1.0)
        parser.add_argument('--length', type=float, default=0.05)
        parser.add_argument('--alpha', type=float, default=1.27e-4)
        
        self.args = parser.parse_args()
        