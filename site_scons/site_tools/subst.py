# File:         subst.py
# Author:       Brian A. Vanderburg II
# Purpose:      A generic SCons file substitution mechanism
# Copyright:    This file is placed in the public domain.
# URL:          http://www.scons.org/wiki/GenericSubstBuilder
##############################################################################


# Requirements
##############################################################################
import re

from SCons.Script import *
import SCons.Errors


# Helper/core functions
##############################################################################

# Do the substitution
def _subst_file(target, source, env, pattern, replace):
    # Read file
    #print 'CALLING SUBST_FILE'
    f = open(source, "rU")
    try:
        contents = f.read()
    finally:
        f.close()

    # Substitute, make sure result is a string
    def subfn(mo):
        value = replace(env, mo)
        if not SCons.Util.is_String(value):
            raise SCons.Errors.UserError("Substitution must be a string.")
        return value
    #print 'pattern  = ' , pattern
    contents = re.sub(pattern, subfn, contents)

    # Write file
    f = open(target, "wt")
    try:
        f.write(contents)
    finally:
        f.close()

# Determine which keys are used
def _subst_keys(source, pattern):
    # Read file
    f = open(source, "rU")
    try:
        contents = f.read()
    finally:
        f.close()

    # Determine keys
    keys = []
    def subfn(mo):
        key = mo.group("key")
        if key:
            keys.append(key)
        return ''

    
    re.sub(pattern, subfn, contents)

    return keys

# Get the value of a key as a string, or None if it is not in the environment
def _subst_value(env, key):
    # Why does "if key in env" result in "KeyError: 0:"?
    try:
        env[key]
    except KeyError:
        return None

    # env.subst already returns a string even if it is stored as a number
    # such as env['HAVE_XYZ'] = 1
    #print 'key = ', key
    #print '  straight env = ', env[key]
    #print '  str of the thing = ', str(env[key])
    #print '  subst(${}) of the thing = ', env.subst("${%s}" % key) 
    #print '  %s of the thing = ', "%s" % str(env[key])
    aa = env[key]
    if aa == []:
       aa = ''
    return aa
    #return str(env[key])
    #return env.subst("${%s}" % key)


# Builder related functions
##############################################################################

# Builder action
def _subst_action(target, source, env):
    # Substitute in the files
    pattern = env["SUBST_PATTERN"]
    replace = env["SUBST_REPLACE"]
    #print 'SUBSTITUTE: ', pattern, ' for ', replace

    for (t, s) in zip(target, source):
        _subst_file(str(t), str(s), env, pattern, replace)

    return 0

# Builder message
def _subst_message(target, source, env):
    items = ["Substituting vars from %s to %s" % (s, t)
             for (t, s) in zip(target, source)]

    return "\n".join(items)

# Builder dependency emitter
def _subst_emitter(target, source, env):
    pattern = env["SUBST_PATTERN"]
    for (t, s) in zip(target, source):
        # When building, if a variant directory is used and source files
        # are being duplicated, the source file will not be duplicated yet
        # when this is called, so the source node must be used instead of
        # the duplicated node
        path = s.srcnode().abspath

        # Get keys used
        keys = _subst_keys(path, pattern)

        d = dict()
        for key in keys:
            value = _subst_value(env, key)
            # print 'key = ', key, ' -> value = ', value
            if not value is None:
                d[key] = value

        # Only the current target depends on this dictionary
        Depends(t, SCons.Node.Python.Value(d))

    return target, source


# Replace @key@ with the value of that key, and @@ with a single @
##############################################################################

_SubstFile_pattern = "@(?P<key>\w*?)@"
def _SubstFile_replace(env, mo):
    key = mo.group("key")
    if not key:
        return "@"

    value = _subst_value(env, key)
    if value is None:
        raise SCons.Errors.UserError("Error: key %s does not exist" % key)
    return value

def SubstFile(env, target, source):
    return env.SubstGeneric(target,
                            source,
                            SUBST_PATTERN=_SubstFile_pattern,
                            SUBST_REPLACE=_SubstFile_replace)


# A substitutor similar to config.h header substitution
# Supported patterns are:
#
# Pattern: #define @key@
# Found:   #define key value
# Missing: /* #define key */
#
# Pattern: #define @key@ default
# Found:   #define key value
# Missing: #define key default
#
# Pattern: #undef @key@
# Found:   #define key value
# Missing: #undef key
#
# The "@" is used to that these defines can be used in addition to
# other defines that you do not desire to be replaced.
##############################################################################

_SubstHeader_pattern = "(?m)^(?P<space>\\s*?)(?P<type>#define|#undef)\\s+?@(?P<key>\w+?)@(?P<ending>.*?)$"
def _SubstHeader_replace(env, mo):
    space = mo.group("space")
    type = mo.group("type")
    key = mo.group("key")
    ending = mo.group("ending")

    value = _subst_value(env, key)
    if not value is None:
        # If found it is always #define key value
        return "%s#define %s %s" % (space, key, value)

    # Not found
    if type == "#define":
        defval = ending.strip()
        if defval:
            # There is a default value
            return "%s#define %s %s" % (space, key, defval)
        else:
            # There is no default value
            return "%s/* #define %s */" % (space, key)

    # It was #undef
    return "%s#undef %s" % (space, key)

def SubstHeader(env, target, source):
    return env.SubstGeneric(target,
                            source,
                            SUBST_PATTERN=_SubstHeader_pattern,
                            SUBST_REPLACE=_SubstHeader_replace)


def generate(env, **kw):
    # The generic builder
    subst = SCons.Action.Action(_subst_action, _subst_message)
    env['BUILDERS']['SubstGeneric'] = Builder(action=subst,
                                              emitter=_subst_emitter)

    # Additional ones
    env.AddMethod(SubstFile)
    env.AddMethod(SubstHeader)


def exists(env):
    return True
