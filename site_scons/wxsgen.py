from __future__ import print_function
import os, sys
import uuid
import xml.etree.ElementTree as et

class WxsGenerator(object):
    def __init__(self, stageDir, short_version, full_version,
                 includeMatlab=False, x64=False):
        self.prefix = stageDir
        self.x64 = x64
        self.includeMatlab = includeMatlab
        self.short_version = short_version
        self.full_version = full_version

        # Use separate UUIDs for 64- and 32-bit components
        if self.x64:
            self.CANTERA_UUID = uuid.UUID('F707EB9E-3723-11E1-A99F-525400631BAF')
            self.pfilesName = 'ProgramFiles64Folder'
            self.productName = 'Cantera {0} (64-bit)'.format(self.short_version)
        else:
            self.CANTERA_UUID = uuid.UUID('1B36CAF0-279D-11E1-8979-001FBC085391')
            self.pfilesName = 'ProgramFilesFolder'
            self.productName = 'Cantera {0} (32-bit)'.format(self.short_version)

    def Directory(self, parent, Id, Name):
        return et.SubElement(parent, 'Directory',
                             dict(Id=Id, Name=Name))

    def FileComponent(self, parent, componentId, fileId, Name, Source,
                      DiskId='1', KeyPath='yes'):
        guid = str(uuid.uuid5(self.CANTERA_UUID, componentId))

        fields = {'Win64': 'yes'} if self.x64 else {}
        c = et.SubElement(parent, "Component",
                          dict(Id=componentId, Guid=guid, **fields))

        fields = {'ProcessorArchitecture': 'x64'} if self.x64 else {}
        f = et.SubElement(c, "File",
                          dict(Id=fileId,
                               Name=Name,
                               Source=Source,
                               DiskId=DiskId,
                               KeyPath=KeyPath,
                               **fields))
        return c,f

    def addDirectoryContents(self, directory, parent, feature):
        """
        directory: name of the directory to add
        parent: the Element for the parent directory
        feature: the Element for the feature to add the files to
        """
        #self.prefix: path to the parent directory
        directories = {}

        # replace characters that are not valid in IDs
        clean = lambda s: s.replace('/', '_').replace('@', 'a').replace('-','_')

        directories[directory] = self.Directory(parent, directory, directory)
        for path, dirs, files in os.walk('/'.join((self.prefix, directory))):
            path = path.replace(self.prefix + '/', '', 1).replace('\\', '/')
            for d in dirs:
                dpath = '/'.join((path, d))
                ID = clean(dpath)
                directories[dpath] = self.Directory(directories[path], ID, d)

            for f in files:
                ID = clean('_'.join((path, f)))
                self.FileComponent(directories[path], ID, ID, f,
                              '/'.join((self.prefix, path, f)))
                et.SubElement(feature, 'ComponentRef', dict(Id=ID))

        return directories

    def make_wxs(self, outFile):
        wix = et.Element("Wix", {'xmlns': 'http://schemas.microsoft.com/wix/2006/wi'})
        product = et.SubElement(wix, "Product",
                                dict(Name=self.productName,
                                     Id=str(self.CANTERA_UUID),
                                     UpgradeCode='2340BEE1-279D-11E1-A4AA-001FBC085391',
                                     Language='1033',
                                     Codepage='1252',
                                     Version=self.full_version,
                                     Manufacturer='Cantera Developers'))

        fields = {'Platform': 'x64'} if self.x64 else {}
        package = et.SubElement(product, "Package",
                                dict(Id='*',
                                     Keywords='Installer',
                                     Description="Cantera {0} Installer".format(self.short_version),
                                     InstallerVersion='310',
                                     Languages='1033',
                                     Compressed='yes',
                                     SummaryCodepage='1252', **fields))

        # Required boilerplate referring to nonexistent installation media
        media = et.SubElement(product, "Media",
                              dict(Id='1',
                                   Cabinet='cantera.cab',
                                   EmbedCab='yes',
                                   DiskPrompt='CD-ROM #1'))
        diskprompt = et.SubElement(product, "Property",
                                   dict(Id='DiskPrompt',
                                        Value="Cantera Installation Disk"))

        # Directories
        targetdir = self.Directory(product, 'TARGETDIR', 'SourceDir')
        pfiles = self.Directory(targetdir, self.pfilesName, 'PFiles')
        instdir = self.Directory(pfiles, 'INSTALLDIR', 'Cantera')

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
        samples = et.SubElement(product, 'Feature',
                               dict(Id='Samples', Level='1',
                                    Title='Samples',
                                    Description='Samples which show you some ways of using Cantera.',
                                    Display='expand',
                                    AllowAdvertise='no'))

        if self.includeMatlab:
            matlab = et.SubElement(product, 'Feature',
                                   dict(Id='Matlab', Level='1',
                                        Title='Matlab Toolbox',
                                        Description='Cantera Matlab toolbox',
                                        Display='expand',
                                        AllowAdvertise='no'))

        # Files
        includes = self.addDirectoryContents('include', instdir, devel)
        binaries = self.addDirectoryContents('bin', instdir, core)
        lib_dir = self.addDirectoryContents('lib', instdir, devel)
        data_dir = self.addDirectoryContents('data', instdir, core)
        sample_dir = self.addDirectoryContents('samples', instdir, samples)
        if self.includeMatlab:
            matlab_dir = self.addDirectoryContents('matlab', instdir, matlab)

        # Registry entries
        reg_options = dict(ForceCreateOnInstall="yes", ForceDeleteOnUninstall="yes",
                           Id='CanteraRegRoot', Root='HKLM',
                           Key='Software\\Cantera\\Cantera {0}'.format(self.short_version))
        reg_key = self.addRegistryKey(core, product, options=reg_options)
        et.SubElement(reg_key, 'RegistryValue', dict(Type='string',
                                                     Name='InstallDir',
                                                     Value='[INSTALLDIR]'))

        # Wix UI
        et.SubElement(product, 'UIRef', dict(Id='WixUI_FeatureTree'))
        et.SubElement(product, 'UIRef', dict(Id='WixUI_ErrorProgressText'))
        et.SubElement(product, 'Property', dict(Id='WIXUI_INSTALLDIR',
                                                Value='INSTALLDIR'))

        # License
        et.SubElement(product, 'WixVariable',
                      dict(Id='WixUILicenseRtf', Value='build/ext/LICENSE.rtf'))

        # Format and save as XML
        indent(wix)
        tree = et.ElementTree(wix)
        tree.write(outFile)

    def addRegistryKey(self, feature, parent, options):
        Id = options['Id']
        guid = str(uuid.uuid5(self.CANTERA_UUID, Id))
        fields = {'Win64': 'yes'} if self.x64 else {}
        dr = et.SubElement(parent, "DirectoryRef", dict(Id="TARGETDIR"))
        c = et.SubElement(dr, "Component", dict(Id=Id, Guid=guid, **fields))
        r = et.SubElement(c, "RegistryKey", options)
        et.SubElement(feature, 'ComponentRef', dict(Id=Id))
        return r


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


def usage():
    print("Usage: wxsgen <stageDir> <outputFile>")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        usage()
        sys.exit()

    WxsGenerator(sys.argv[1]).make_wxs(sys.argv[2])
