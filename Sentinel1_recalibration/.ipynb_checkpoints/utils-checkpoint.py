import os
import shutil



def get_memory_usage(unit='Go', var='ru_maxrss', force_psutil=False):
    """
    var str: ru_maxrss or ru_ixrss or ru_idrss or ru_isrss or current
    Returns
    -------

    """
    if unit == 'Go':
        factor = 1000000.
    elif unit == 'Mo':
        factor = 1000.
    elif unit == 'Ko':
        factor = 1.
    else:
        raise Exception('not handle unit')

    try:
        if force_psutil:
            # on_purpose_error
            None
        import resource
        mems = {}
        mems['ru_maxrss'] = resource.getrusage(
            resource.RUSAGE_SELF).ru_maxrss / factor
        mems['ru_ixrss'] = resource.getrusage(
            resource.RUSAGE_SELF).ru_ixrss / factor
        mems['ru_idrss'] = resource.getrusage(
            resource.RUSAGE_SELF).ru_idrss / factor
        mems['ru_isrss'] = resource.getrusage(
            resource.RUSAGE_SELF).ru_isrss / factor
        mems['current'] = getCurrentMemoryUsage() / factor
        # memory_used_go = resource.getrusage(resource.RUSAGE_SELF).get(var) /factor
        memory_used_go = mems[var]
    except:  # on windows resource is not usable
        import psutil
        memory_used_go = psutil.virtual_memory().used / factor / 1000.
    str_mem = 'RAM usage: %1.1f %s' % (memory_used_go, unit)
    return str_mem


def getCurrentMemoryUsage():
    ''' Memory usage in kB '''

    with open('/proc/self/status') as f:
        memusage = f.read().split('VmRSS:')[1].split('\n')[0][:-3]

    return int(memusage.strip())


def copy_tree(src, dst, ignore=None):
    if not os.path.exists(dst):
        os.makedirs(dst)
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            copy_tree(s, d, ignore)
        else:
            if not ignore or not ignore(s):
                shutil.copy2(s, d)
                
def ignore_files(file):
    #ignore tiff
    return file.endswith('.tiff')