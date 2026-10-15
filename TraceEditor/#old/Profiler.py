import cProfile, pstats

def profile(func):
    '''
    Decorator for profiling a function.

    Parameters
    ----------
    func : function
        The function to be profiled.

    '''
    def wrapper(*args, **kwargs):
        filename = f'''profile_{func.__name__}.pstats'''
        print(f'>>>>>> Profiling {func.__name__}')
        profiler = cProfile.Profile()
        profiler.enable()
        
        result = func(*args, **kwargs)
        
        profiler.disable()
        # profiler.dump_stats(filename)
        print(f'>>>>>> Profiling results: {filename}')
        stats = pstats.Stats(filename)
        stats.sort_stats('cumulative').print_stats(10) 

        return result
    
    return wrapper