from math          import log, ceil
from scipy.special import betaincinv as ribf

import random

zero_const = 1. / 10**10

def get_trials(prob):
    return 2 * int(ceil(log(zero_const) / log(1. - prob)))

def get_prob(trials, success):
    # probability for rare event success
    return (success + 1.) / (trials + 2.)

def get_CI(trials, success, conf_interval):
    x1 = ribf(success + 1, trials - success + 1, 0.5*(1 - conf_interval))
    x2 = ribf(success + 1, trials - success + 1, 0.5*(1 + conf_interval))
    return (x1, x2)

def test_func():
    if random.random() <= 0.002: #1. / 200:
        return True
    return False

def main():
    curr_prob     = 1. / 100000
    conf_interval = 0.9999

    total_success = 1.
    total_trials  = get_init_trials(recommended_trials=50)
    curr_prob     = get_prob(total_trials, total_success)
    new_trial     = total_trials
    
    attempt       = 0

    while True:

        for t in xrange(new_trial):
            total_trials += 1
            if test_func():
                total_success += 1.
                print '[{}] Current Prob = {}'.format(attempt, get_prob(total_trials, total_success))
                attempt += 1
                break
        else:
            print
            print '  Est Prob = {}'.format(get_prob(total_trials, total_success))
            x1, x2 = get_CI(total_trials, total_success, conf_interval)
            print '-99.99% CI = {}'.format(x1)
            print '+99.99% CI = {}'.format(x2)
            return None

        curr_prob = get_prob(total_trials, total_success)
        new_trial = get_trials(curr_prob)

if __name__ == '__main__':
    main()