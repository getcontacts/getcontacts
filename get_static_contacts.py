#!/usr/bin/env python3

import dynamic_contacts
import os

if __name__ == "__main__":
    dynamic_contacts.main(traj_required=False)

    # Suppress stdout from vmd as program terminates
    devnull = open('/dev/null', "w")
    os.dup2(devnull.fileno(), 1)
