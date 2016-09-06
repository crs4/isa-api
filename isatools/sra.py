import logging
import os
import iso8601
import jinja2
import html
import datetime
import hashlib
from functools import partial
from lxml import etree
import xml.dom.minidom
from isatools.model.v1 import Sample, OntologyAnnotation, DataFile

logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)

supported_sra_assays = [
    ('genome sequencing', 'nucleotide sequencing'),
    ('environmental gene survey', 'nucleotide sequencing'),
    ('metagenome sequencing', 'nucleotide sequencing')
]

sra_center_name = 'OXFORD'
sra_broker_name = 'ISAcreator'
sra_lab_name = sra_center_name
sra_submission_action = 'ADD'
sra_center_prj_name = None


def export(investigation, export_path, sra_settings=None, datafilehashes=None):

    def get_comment(assay, name):
        hits = [c for c in assay.comments if c.name.lower() == name.lower()]
        if len(hits) > 1:
            raise AttributeError("Multiple comments of label '{}' found".format(name))
        elif len(hits) < 1:
            return None
        else:
            return hits[0]

    def get_sample(process):
        materials = process.inputs
        sample = None
        for material in materials:
            if isinstance(material, Sample):
                sample = material
                break
        return sample

    def get_pv(process, name):
        hits = [pv for pv in process.parameter_values if
                pv.category.parameter_name.name.lower().replace('_', ' ') == name.lower().replace('_', ' ')]
        if len(hits) > 1:
            raise AttributeError("Multiple parameter values of category '{}' found".format(name))
        elif len(hits) < 1:
            return None
        else:
            if isinstance(hits[0].value, OntologyAnnotation):
                value = hits[0].value.name
            else:
                value = hits[0].value
            return value.replace('_', ' ')

    global sra_center_name
    global sra_broker_name
    if sra_settings is not None:
        sra_center_name = sra_settings['sra_center']
        sra_broker_name = sra_settings['sra_broker']
        # sra_lab_name = sra_settings['sra_lab_name']
        # sra_submission_action = sra_settings['sra_submission_action']
        # sra_center_prj_name = sra_settings['sra_center_prj_name']

    logger.info("isatools.sra.export()")
    for istudy in investigation.studies:
        is_sra = False
        for iassay in istudy.assays:
            if (iassay.measurement_type.name, iassay.technology_type.name) in supported_sra_assays:
                is_sra = True
                break
        if not is_sra:
            logger.info("No SRA assay found, skipping processing")
            continue

        study_acc = istudy.identifier
        logger.debug("sra exporter, working on " + study_acc)

        # Flag SRA contacts for template
        has_sra_contact = False
        for contact in istudy.contacts:
            if "sra inform on status" in [r.name.lower() for r in contact.roles]:
                contact.inform_on_status = True
                has_sra_contact = True
            if "sra inform on error" in [r.name.lower() for r in contact.roles]:
                contact.inform_on_error = True
                has_sra_contact = True
        if not has_sra_contact:
            raise ValueError(
                "The study ''{0}'' has either no SRA contact or no email specified for the contact. Please "
                "ensure you have one contact with a 'Role' as 'SRA Inform On Status', otherwise we cannot "
                "export to SRA.".format(istudy.identifier))

        if istudy.submission_date is None or istudy.submission_date == '':
            istudy.submission_date = iso8601.parse_date(datetime.date.today().isoformat(), iso8601.UTC).isoformat()
        else:
            istudy.submission_date = iso8601.parse_date(istudy.submission_date, iso8601.UTC).isoformat()
        istudy.description = html.escape(istudy.description)  # ideally make it a requirement in the model or JSON to have html escaped content

        env = jinja2.Environment()
        env.loader = jinja2.FileSystemLoader(os.path.join(os.path.dirname(__file__), 'resources', 'sra_templates'))
        xsub_template = env.get_template('submission_add.xml')
        xsub = xsub_template.render(accession=study_acc, contacts=istudy.contacts, submission_date=istudy.submission_date,
                                    sra_center_name=sra_center_name, sra_broker_name=sra_broker_name)
        xproj_template = env.get_template('project_set.xml')
        xproj = xproj_template.render(study=istudy, sra_center_name=sra_center_name)

        assays_to_export = list()
        for iassay in istudy.assays:
            if (iassay.measurement_type.name, iassay.technology_type.name) in supported_sra_assays:
                assay_seq_processes = [a for a in iassay.process_sequence if a.executes_protocol.protocol_type.name ==
                                       'nucleic acid sequencing']
                for assay_seq_process in assay_seq_processes:
                    do_export = True
                    if get_comment(assay_seq_process, 'export') is not None:
                        logger.debug("HAS EXPORT COMMENT IN ASSAY")
                        export = get_comment(assay_seq_process, 'export').value
                        logger.debug("export is " + export)
                        do_export = export.lower() != 'no'
                    else:
                        logger.debug("NO EXPORT COMMENT FOUND")
                    logger.debug("Perform export? " + str(do_export))
                    if do_export:
                        sample = None
                        curr_process = assay_seq_process
                        while sample is None:
                            sample = get_sample(curr_process)
                            curr_process = curr_process.prev_process
                        assay_to_export = \
                            {
                                "sample": sample,
                                "sample_alias": study_acc + ':sample:' + sample.name,
                                "run_alias": study_acc + ':assay:' + assay_seq_process.name,
                                "exp_alias": study_acc + ':generic_assay:' + iassay.filename[:-4] + ':' + assay_seq_process.name
                            }
                        datafiles = [d for d in assay_seq_process.outputs if isinstance(d, DataFile)]
                        checksum = '00000000000000000000000000000000'
                        if datafilehashes is not None:
                            checksum = datafilehashes[datafiles[0].filename]  # raises AttributeError if file not found
                        filetype = datafiles[0].filename[datafiles[0].filename.index('.') + 1:]
                        if filetype.endswith('.gz'):
                            filetype = filetype[:filetype.index('.')]
                        assay_to_export['data_file'] = {
                            "filename": datafiles[0].filename,
                            "filetype": filetype,
                            "checksum": checksum
                        }
                        source = None
                        matching_sources = [p.inputs for p in istudy.process_sequence if sample in p.outputs]
                        if len(matching_sources[0]) == 1:
                            source = matching_sources[0][0]
                        assay_to_export['source'] = {
                            "name": source.name,
                            "characteristics": source.characteristics,
                        }
                        organism_charac = [c for c in source.characteristics if c.category.name == 'organism'][-1]
                        assay_to_export['source']['taxon_id'] = organism_charac.value.term_accession[organism_charac.value.term_accession.index('_')+1:]
                        assay_to_export['source']['scientific_name'] = organism_charac.value.name
                        curr_process = assay_seq_process
                        while curr_process.prev_process is not None:
                            assay_to_export[curr_process.executes_protocol.protocol_type.name] = curr_process
                            try:
                                curr_process = curr_process.prev_process
                            except AttributeError:
                                pass
                        target_taxon = get_pv(assay_to_export['library construction'], 'target_taxon')
                        assay_to_export['target_taxon'] = target_taxon
                        assay_to_export['targeted_loci'] = False
                        # BEGIN genome seq library selection
                        if iassay.measurement_type.name in ['genome sequencing', 'whole genome sequencing']:
                            library_source = get_pv(assay_to_export['library construction'],
                                                      'library source')
                            if library_source.upper() not in ['GENOMIC', 'GENOMIC SINGLE CELL', 'METAGENOMIC', 'OTHER']:
                                logger.warn("ERROR:value supplied is not compatible with SRA1.5 schema " + library_source)
                                library_source = 'OTHER'

                            library_strategy = get_pv(assay_to_export['library construction'],
                                                      'library strategy')
                            if library_strategy.upper() not in ['WGS', 'OTHER']:
                                logger.warn("ERROR:value supplied is not compatible with SRA1.5 schema " + library_strategy)
                                library_strategy = 'OTHER'

                            library_selection = get_pv(assay_to_export['library construction'],
                                                       'library selection')
                            if library_selection not in ['RANDOM', 'UNSPECIFIED']:
                                logger.warn("ERROR:value supplied is not compatible with SRA1.5 schema " + library_selection)
                                library_selection = 'unspecified'

                            protocol = "\n protocol_description: " \
                                       + assay_to_export['library construction'].executes_protocol.description
                            mid_pv = get_pv(assay_to_export['library construction'], 'mid')
                            if mid_pv is not None:
                                protocol += "\n mid: " + mid_pv.value

                            assay_to_export['library_source'] = library_source
                            assay_to_export['library_strategy'] = library_strategy
                            assay_to_export['library_selection'] = library_selection
                            assay_to_export['library_construction_protocol'] = protocol

                            library_layout = get_pv(assay_to_export['library construction'], 'library layout')
                            assay_to_export['library_layout'] = library_layout
                        # END genome seq library selection
                        # BEGIN environmental gene survey library selection
                        elif iassay.measurement_type.name in ['environmental gene survey']:
                            assay_to_export['library_source'] = 'METAGENOMIC'
                            assay_to_export['library_strategy'] = 'AMPLICON'
                            assay_to_export['library_selection'] = 'PCR'
                            library_layout = get_pv(assay_to_export['library construction'], 'library layout')
                            assay_to_export['library_layout'] = library_layout
                            nucl_acid_amp = get_pv(assay_to_export['library construction'], 'nucleic acid amplification')
                            if nucl_acid_amp is None:
                                nucl_acid_amp = get_pv(assay_to_export['library construction'], 'nucl_acid_amp')

                            protocol = "\n protocol_description: " \
                                       + assay_to_export['library construction'].executes_protocol.description
                            mid_pv = get_pv(assay_to_export['library construction'], 'mid')
                            if mid_pv is not None:
                                protocol += "\n mid: " + mid_pv
                                assay_to_export['barcode'] = mid_pv
                                assay_to_export['min_match'] = len(mid_pv)
                            if nucl_acid_amp is not None:
                                protocol += "\n nucl_acid_amp: " + nucl_acid_amp.value
                            url = get_pv(assay_to_export['library construction'], 'url')
                            if url is not None:
                                protocol += "\n url: " + nucl_acid_amp.value
                            target_taxon = assay_to_export['target_taxon']
                            if target_taxon is not None:
                                protocol += "\n target_taxon: " + target_taxon
                            target_gene = get_pv(assay_to_export['library construction'], 'target_gene')
                            if target_gene is not None:
                                protocol += "\n target_gene: " + target_gene
                            target_subfragment = get_pv(assay_to_export['library construction'], 'target_subfragment')
                            if target_subfragment is not None:
                                protocol += "\n target_subfragment: " + target_subfragment
                            pcr_primers = get_pv(assay_to_export['library construction'], 'pcr_primers')
                            if pcr_primers is not None:
                                protocol += "\n pcr_primers: " + pcr_primers
                            pcr_cond = get_pv(assay_to_export['library construction'], 'pcr_cond')
                            if pcr_cond is not None:
                                protocol += "\n pcr_cond: " + pcr_cond
                            assay_to_export['library_construction_protocol'] = protocol

                            if target_gene is not None:
                                assay_to_export['targeted_loci'] = True
                                assay_to_export['locus_name'] = target_gene
                        # END environmental gene survey library selection
                        # BEGIN metagenome seq library selection
                        elif iassay.measurement_type.name in ['metagenome sequencing']:
                            library_source = 'METAGENOMIC'
                            library_strategy = get_pv(assay_to_export['library construction'],
                                                      'library strategy')
                            if library_strategy.upper() not in ['WGS', 'OTHER']:
                                logger.warn(
                                    "ERROR:value supplied is not compatible with SRA1.5 schema " + library_strategy)
                                library_strategy = 'OTHER'

                            library_selection = get_pv(assay_to_export['library construction'],
                                                       'library selection')
                            if library_selection not in ['RANDOM', 'UNSPECIFIED']:
                                logger.warn(
                                    "ERROR:value supplied is not compatible with SRA1.5 schema " + library_selection)
                                library_selection = 'unspecified'

                            protocol = "\n protocol_description: " \
                                       + assay_to_export['library construction'].executes_protocol.description
                            mid_pv = get_pv(assay_to_export['library construction'], 'mid')
                            if mid_pv is not None:
                                protocol += "\n mid: " + mid_pv.value

                            assay_to_export['library_source'] = library_source
                            assay_to_export['library_strategy'] = library_strategy
                            assay_to_export['library_selection'] = library_selection
                            assay_to_export['library_construction_protocol'] = protocol

                            library_layout = get_pv(assay_to_export['library construction'], 'library layout')
                            assay_to_export['library_layout'] = library_layout
                        # END metagenome seq library selection
                        # BEGIN transciption profiling library selection
                        elif iassay.measurement_type.name in ['transcription profiling']:
                            # TODO: Implement this section to get BII-S-3 working
                            pass
                        # END transciption profiling library selection
                        else:
                            logger.error("ERROR:Unsupported measurement type: " + iassay.measurement_type.name)
                        assay_to_export['platform'] = get_pv(assay_to_export['nucleic acid sequencing'],
                                                             'sequencing instrument')
                        assays_to_export.append(assay_to_export)
            else:
                logger.error("ERROR:Unsupported measurement/technology type {0}/{1}, skipping assays".format(iassay.measurement_type.name, iassay.technology_type.name))

        xexp_set_template = env.get_template('experiment_set.xml')
        xexp_set = xexp_set_template.render(assays_to_export=assays_to_export, study=istudy,
                                            sra_center_name=sra_center_name, sra_broker_name=sra_broker_name)
        xrun_set_template = env.get_template('run_set.xml')
        xrun_set = xrun_set_template.render(assays_to_export=assays_to_export, study=istudy,
                                            sra_center_name=sra_center_name, sra_broker_name=sra_broker_name)

        xsample_set_template = env.get_template('sample_set.xml')
        xsample_set = xsample_set_template.render(assays_to_export=assays_to_export, study=istudy,
                                            sra_center_name=sra_center_name, sra_broker_name=sra_broker_name)
        logger.debug("SRA exporter: writing SRA XML files for study " + study_acc)

        # blitz out whitespaces with etree and format nicely with minidom
        def prettify(xmlstr):
            p = etree.XMLParser(remove_blank_text=True)
            exsub = etree.XML(xmlstr, parser=p)
            x = xml.dom.minidom.parseString(etree.tostring(exsub))
            return x.toprettyxml()

        # validate a doc against one of our SRA schemas
        def validate(docpath, schemaname):
            with open(os.path.join(os.path.dirname(__file__), 'resources', 'sra_schemas', schemaname)) as xsd:
                doc = etree.parse(xsd)
                try:
                    schema = etree.XMLSchema(doc)
                    with open(docpath, 'r') as xsub_file:
                        doc = etree.parse(xsub_file)
                        try:
                            schema.assertValid(doc)
                        except etree.DocumentInvalid as e:
                            logger.error("Schema validation failed on " + docpath + ':\n' + str(e))
                except etree.XMLSchemaParseError as e:
                    logger.error(e)

        if os.path.exists(export_path):
            with open(os.path.join(export_path, 'submission.xml'), 'w') as xsub_file:
                print(prettify(xsub), file=xsub_file)
            validate(os.path.join(export_path, 'submission.xml'), 'SRA.submission.xsd')
            with open(os.path.join(export_path, 'project_set.xml'), 'w') as xproj_set_file:
                print(prettify(xproj), file=xproj_set_file)
            validate(os.path.join(export_path, 'project_set.xml'), 'ENA.project.xsd')
            with open(os.path.join(export_path, 'experiment_set.xml'), 'w') as xexp_set_file:
                print(prettify(xexp_set), file=xexp_set_file)
            validate(os.path.join(export_path, 'experiment_set.xml'), 'SRA.experiment.xsd')
            with open(os.path.join(export_path, 'run_set.xml'), 'w') as xrun_set_file:
                print(prettify(xrun_set), file=xrun_set_file)
            validate(os.path.join(export_path, 'run_set.xml'), 'SRA.run.xsd')
            with open(os.path.join(export_path, 'sample_set.xml'), 'w') as xsample_set_file:
                print(prettify(xsample_set), file=xsample_set_file)
            validate(os.path.join(export_path, 'sample_set.xml'), 'SRA.sample.xsd')
        else:
            raise NotADirectoryError("export path '{}' is not a directory".format(export_path))


def create_datafile_hashes(fileroot, filenames):
    """
    Create md5 file dict for files in a directory with a particular extension

    :param fileroot: Root to directory containing files (assumes all in same dir)
    :param filenames: List of filenames of files to md5, assumed in fileroot
    :return: dict containing filenames and md5s

    Usage:
    >>> filenames = [f for f in listdir('/path/to/my/files') if f.endswith('.fastq.gz')]
    >>> create_datafile_hashes(fileroot='/path/to/my/files', filenames=filesnames)
    { 'myfile1.gz': 'd41d8cd98f00b204e9800998ecf8427e', 'myfile2.gz': 'd41d8cd98f00b204e9800998ecf8427e' }
    """
    def md5sum(filename):
        with open(filename, mode='rb') as f:
            d = hashlib.md5()
            for buf in iter(partial(f.read, 128), b''):
                d.update(buf)
        return d.hexdigest()

    from os.path import isfile, join
    datafilehashes = dict()
    for file in filenames:
        if isfile(join(fileroot, file)):
            datafilehashes[file] = md5sum(filename=join(fileroot, file))
        else:
            raise FileNotFoundError("{} is not a file".format(fileroot + file))
    return datafilehashes