try:
    # Prefer LXML for its pretty printer, but xml.etree works fine.
    import lxml.etree as et
    have_lxml = True
except ImportError:
    import xml.etree.ElementTree as et
    have_lxml = False

wix = et.Element("Wix", {'xmlns': 'http://schemas.microsoft.com/wix/2006/wi'})

def Directory(parent, Id, Name):
    return et.SubElement(parent, 'Directory',
                         dict(Id=Id, Name=Name))

def FileComponent(parent, componentId, fileId, Guid, Name, Source, DiskId='1', KeyPath='yes'):
    c = et.SubElement(parent, "Component",
                      dict(Id=componentId, Guid=Guid))
    f = et.SubElement(c, "File",
                      dict(Id=fileId,
                           Name=Name,
                           Source=Source,
                           DiskId=DiskId,
                           KeyPath=KeyPath))
    return c,f

product = et.SubElement(wix, "Product",
                     dict(Name='Cantera 2.0',
                          Id='1B36CAF0-279D-11E1-8979-001FBC085391',
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

targetdir = Directory(product, 'TARGETDIR', 'SourceDir')
pfiles = Directory(targetdir, 'ProgramFilesFolder', 'PFiles')
instdir = Directory(pfiles, 'INSTALLDIR', 'Cantera')
bindir = Directory(instdir, 'bin', 'bin')

ck2cti = FileComponent(bindir, 'MainExecutable', 'ck2ctiEXE', 
                       '0ABCA730-279D-11E1-984E-001FBC085391',
                       'ck2cti.exe', '../../stage/bin/ck2cti.exe')

complete = et.SubElement(product, 'Feature', dict(Id='Complete', Level='1'))
mainExec = et.SubElement(complete, 'ComponentRef', dict(Id='MainExecutable'))

if have_lxml:
    print et.tostring(wix, pretty_print=True)
else:
    print et.tostring(wix)
