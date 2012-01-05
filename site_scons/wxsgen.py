import os, sys
import uuid
import xml.etree.ElementTree as et

CANTERA_UUID = uuid.UUID('1B36CAF0-279D-11E1-8979-001FBC085391')

def Directory(parent, Id, Name):
    return et.SubElement(parent, 'Directory',
                         dict(Id=Id, Name=Name))

def FileComponent(parent, componentId, fileId, Name, Source, DiskId='1', KeyPath='yes'):
    guid = str(uuid.uuid5(CANTERA_UUID, componentId))
    c = et.SubElement(parent, "Component",
                      dict(Id=componentId, Guid=guid))
    f = et.SubElement(c, "File",
                      dict(Id=fileId,
                           Name=Name,
                           Source=Source,
                           DiskId=DiskId,
                           KeyPath=KeyPath))
    return c,f


def indent(elem, level=0):
    """ in-place prettyprint formatter (from lxml) """
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


def addDirectoryContents(prefix, directory, parent, feature):
    """
    prefix: path to the parent directory
    directory: name of the directory to add
    parent: the Element for the parent directory
    feature: the Element for the feature to add the files to
    """
    directories = {}

    directories[directory] = Directory(parent, directory, directory)
    for path, dirs, files in os.walk('/'.join((prefix, directory))):
        path = path.replace(prefix + '/', '', 1).replace('\\', '/')
        for d in dirs:
            dpath = '/'.join((path, d))
            ID = dpath.replace('/', '_')
            directories[dpath] = Directory(directories[path], ID, d)

        for f in files:
            ID = '_'.join((path, f)).replace('/', '_')
            FileComponent(directories[path], ID, ID, f,
                          '/'.join((prefix, path, f)))
            et.SubElement(feature, 'ComponentRef', dict(Id=ID))

    return directories


def make_wxs(stageDir, outFile):
    wix = et.Element("Wix", {'xmlns': 'http://schemas.microsoft.com/wix/2006/wi'})
    product = et.SubElement(wix, "Product",
                            dict(Name='Cantera 2.0',
                                 Id=str(CANTERA_UUID),
                                 UpgradeCode='2340BEE1-279D-11E1-A4AA-001FBC085391',
                                 Language='1033',
                                 Codepage='1252',
                                 Version='2.0.0',
                                 Manufacturer='Cantera Developers'))

    package = et.SubElement(product, "Package",
                            dict(Id='*',
                                 Keywords='Installer',
                                 Description="Cantera 2.0 Installer",
                                 InstallerVersion='100',
                                 Languages='1033',
                                 Compressed='yes',
                                 SummaryCodepage='1252'))

    # Required boilerplate refering to nonexistent installation media
    media = et.SubElement(product, "Media",
                          dict(Id='1',
                               Cabinet='cantera.cab',
                               EmbedCab='yes',
                               DiskPrompt='CD-ROM #1'))
    diskprompt = et.SubElement(product, "Property",
                               dict(Id='DiskPrompt',
                                    Value="Cantera Installation Disk"))

    # Directories
    targetdir = Directory(product, 'TARGETDIR', 'SourceDir')
    pfiles = Directory(targetdir, 'ProgramFilesFolder', 'PFiles')
    instdir = Directory(pfiles, 'INSTALLDIR', 'Cantera')

    # Features
    core = et.SubElement(product, 'Feature',
                             dict(Id='Core', Level='1',
                                  Title='Cantera',
                                  Description='Cantera base files',
                                  Display='expand',
                                  ConfigurableDirectory='INSTALLDIR',
                                  AllowAdvertise='no',
                                  Absent='disallow'))
    devel = et.SubElement(product, 'Feature',
                          dict(Id='DevTools', Level='1000',
                               Title='Develpment Tools',
                               Description='Header files and static libraries needed to develop applications that use Cantera.',
                               Display='expand',
                               AllowAdvertise='no'))
    extras = et.SubElement(product, 'Feature',
                           dict(Id='Extras', Level='1',
                                Title='Extras',
                                Description='Demos, tutorials and templates which show you some ways of using Cantera.',
                                Display='expand',
                                AllowAdvertise='no'))

    # Files
    includes = addDirectoryContents(stageDir, 'include', instdir, devel)
    binaries = addDirectoryContents(stageDir, 'bin', instdir, core)
    lib_dir = addDirectoryContents(stageDir, 'lib', instdir, devel)
    data_dir = addDirectoryContents(stageDir, 'data', instdir, core)
    demos_dir = addDirectoryContents(stageDir, 'demos', instdir, extras)
    templates_dir = addDirectoryContents(stageDir, 'templates', instdir, extras)
    tutorials_dir = addDirectoryContents(stageDir, 'tutorials', instdir, extras)

    # Wix UI
    et.SubElement(product, 'UIRef', dict(Id='WixUI_FeatureTree'))
    et.SubElement(product, 'UIRef', dict(Id='WixUI_ErrorProgressText'))
    et.SubElement(product, 'Property', dict(Id='WIXUI_INSTALLDIR',
                                            Value='INSTALLDIR'))

    # License
    et.SubElement(product, 'WixVariable',
                  dict(Id='WixUILicenseRtf', Value='platform/windows/License.rtf'))

    # Format and save as XML
    indent(wix)
    tree = et.ElementTree(wix)
    tree.write(outFile)

def usage():
    print "Usage: wxsgen <stageDir> <outputFile>"

if __name__ == '__main__':
    if len(sys.argv) != 3:
        usage()
        sys.exit()

    make_wxs(sys.argv[1], sys.argv[2])
