class DefineDict(object):
    def __init__(self, data):
        self.data = data
        self.undefined = set()

    def __getitem__(self, key):
        if key not in self.data:
            self.undefined.add(key)
            return '/* #undef %s */' % key
        elif self.data[key] is None:
            return '/* #undef %s */' % key
        else:
            return '#define %s %s' % (key, self.data[key])

class ConfigBuilder(object):
    def __init__(self, defines):
        self.defines = DefineDict(defines)

    def __call__(self, source, target, env):
        for s, t in zip(source, target):
            config_h_in = file(str(s), "r")
            config_h = file(str(t), "w")

            config_h.write(config_h_in.read() % self.defines)
            config_h_in.close()
            config_h.close()
            self.print_config(str(t))

    def print_config(self, filename):
        print 'Generating %s with the following settings:' % filename
        for key, val in sorted(self.defines.data.iteritems()):
            if val is not None:
                print "    %-35s %s" % (key, val)
        for key in sorted(self.defines.undefined):
            print "    %-35s %s" % (key, '*undefined*')

def quoted(s):
    return '"%s"' % s
