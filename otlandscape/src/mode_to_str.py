def mode_to_str(mode):
    if( mode == 1 ):
        return 'W2 trace split-renormalize surface'
    elif( mode == 2 ):
        return 'L2 surface'
    elif( mode == 3 ):
        return 'W2 trace split-renormalize global'
    elif( mode == 4 ):
        return 'L2 global'
    elif( mode == 5 ):
        return 'W2 trace abs-renormalize surface'
    else:
        print('Mode not supported...exiting')
        exit(-1)

