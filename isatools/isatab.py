"""Functions for reading, writing and validating ISA-Tab.

This version only implements
    - ISA-Tab parsing (no graph support yet)
    - Simple investigation and study sample table file validation
    - A basic file sniffer

Don't forget to read the ISA-Tab spec:
http://isa-specs.readthedocs.io/en/latest/isatab.html
"""
from __future__ import absolute_import
from bisect import bisect_left, bisect_right
import csv
import itertools
import logging
import os
import re
from collections import namedtuple
from collections import OrderedDict
from six.moves import zip_longest, zip

from isatools.model import *

__author__ = 'djcomlab@gmail.com (David Johnson)'

logging.basicConfig(level=logging.ERROR)
log = logging.getLogger(__name__)


class Sniffer(object):

    def __init__(self):
        self.is_tab_delimited = None
        self.is_isatab_investigation_filename = None
        self.is_isatab_study_filename = None
        self.is_isatab_assay_filename = None

    def sniff(self, filepath_or_buffer):
        if isinstance(filepath_or_buffer, str):
            with open(filepath_or_buffer, 'rU') as filebuffer:
                return self._sniff(filebuffer)
        else:
            return self._sniff(filepath_or_buffer)

    def _sniff(self, filebuffer):
        filename = os.path.basename(filebuffer.name)
        sniff = csv.Sniffer().sniff(filebuffer.read(1024))
        if sniff.delimiter == '\t':
            self.is_tab_delimited = True
        filebuffer.seek(0)
        if filename.startswith('i_') and filename.endswith('.txt'):
            self.is_isatab_investigation_filename = True
            self.is_isatab_study_filename = False
            self.is_isatab_assay_filename = False
        elif filename.startswith('s_') and filename.endswith('.txt'):
            self.is_isatab_investigation_filename = False
            self.is_isatab_study_filename = True
            self.is_isatab_assay_filename = False
        elif filename.startswith('a_') and filename.endswith('.txt'):
            self.is_isatab_investigation_filename = False
            self.is_isatab_study_filename = False
            self.is_isatab_assay_filename = True


class AbstractParser(object):

    @staticmethod
    def _pairwise(iterable):
        a, b = itertools.tee(iterable)
        next(b, None)
        return zip(a, b)

    def parse(self, filepath_or_buffer):
        if isinstance(filepath_or_buffer, str):
            with open(filepath_or_buffer, 'rU') as filebuffer:
                self._parse(filebuffer)
        else:
            self._parse(filepath_or_buffer)

    def _parse(self, filebuffer):
        raise NotImplementedError(
            'Inherit from this class and implement this method')


class TableParser(AbstractParser):

    DATA_FILE_LABELS = (
        'Raw Data File', 'Derived Spectral Data File',
        'Derived Array Data File', 'Array Data File',
        'Protein Assignment File', 'Peptide Assignment File',
        'Post Translational Modification Assignment File',
        'Acquisition Parameter Data File', 'Free Induction Decay Data File',
        'Derived Array Data Matrix File', 'Image File', 'Derived Data File',
        'Metabolite Assignment File', 'Raw Spectral Data File')
    MATERIAL_LABELS = ('Source Name', 'Sample Name', 'Extract Name',
                       'Labeled Extract Name')
    OTHER_MATERIAL_LABELS = ('Extract Name', 'Labeled Extract Name')
    NODE_LABELS = DATA_FILE_LABELS + MATERIAL_LABELS + OTHER_MATERIAL_LABELS
    ASSAY_LABELS = ('Assay Name', 'MS Assay Name', 'Hybridization Assay Name',
                    'Scan Name', 'Data Transformation Name',
                    'Normalization Name')
    ALL_LABELS = NODE_LABELS + ASSAY_LABELS + tuple('Protocol REF')

    @staticmethod
    def _find_lt(a, x):
        i = bisect_left(a, x)
        if i:
            return a[i - 1]
        else:
            return -1

    @staticmethod
    def _find_gt(a, x):
        i = bisect_right(a, x)
        if i != len(a):
            return a[i]
        else:
            return -1

    @staticmethod
    def _clean_label(label):
        for clean_label in TableParser.ALL_LABELS:
            if label.startswith(clean_label):
                return clean_label

    def __init__(self):
        self.node_map = OrderedDict()
        self.process_map = OrderedDict()

    def _make_process_sequence(self, df):
        df = df[[x for x in df.columns if
                 x.startswith(TableParser.ALL_LABELS)]]
        process_key_sequences = []
        for rindex, row in df.iterrows():
            process_key_sequence = []
            labels = df.columns
            nodes_index = [i for i, x in enumerate(labels) if
                           x in TableParser.NODE_LABELS]
            for cindex, label in enumerate(labels):
                val = row[label]
                if label.startswith('Protocol REF') and val != '':
                    output_node_index = self._find_gt(nodes_index, cindex)
                    if output_node_index > -1:
                        output_node_label = labels[output_node_index]
                        output_node_val = row[output_node_label]
                    input_node_index = self._find_lt(nodes_index, cindex)
                    if input_node_index > -1:
                        input_node_label = labels[input_node_index]
                        input_node_val = row[input_node_label]
                    input_nodes_with_prot_keys = df.loc[
                        df[labels[cindex]] == val].groupby(
                        [labels[cindex], labels[input_node_index]]).size()
                    output_nodes_with_prot_keys = df.loc[
                        df[labels[cindex]] == val].groupby(
                        [labels[cindex], labels[output_node_index]]).size()
                    if len(input_nodes_with_prot_keys) > len(
                            output_nodes_with_prot_keys):
                        process_key = '.'.join([val, output_node_val.strip()])
                    elif len(input_nodes_with_prot_keys) < len(
                            output_nodes_with_prot_keys):
                        process_key = '.'.join([input_node_val.strip(), val])
                    else:
                        process_key = '.'.join([input_node_val.strip(), val,
                                                output_node_val.strip()])
                    if process_key not in self.process_map.keys():
                        process = Process(id_=process_key)
                        self.process_map[process_key] = process
                    process_key_sequence.append(process_key)
                elif label.startswith(TableParser.NODE_LABELS):
                    process_key_sequence.append(
                        '.'.join([self._clean_label(label), val]))
            process_key_sequences.append(process_key_sequence)
        for process_key_sequence in process_key_sequences:
            for left, right in self._pairwise(process_key_sequence):
                if left.startswith(TableParser.NODE_LABELS) and not \
                        right.startswith(TableParser.NODE_LABELS):
                    material = self.node_map[left]
                    process = self.process_map[right]
                    if material not in process.inputs:
                        process.inputs.append(material)
                elif not left.startswith(TableParser.NODE_LABELS) and \
                        right.startswith(TableParser.NODE_LABELS):
                    process = self.process_map[left]
                    material = self.node_map[right]
                    if material not in process.outputs:
                        process.outputs.append(material)
        for process_key_sequence in process_key_sequences:
            process_only_key_sequence = filter(
                lambda x: not x.startswith(TableParser.NODE_LABELS),
                process_key_sequence)
            for left, right in self._pairwise(process_only_key_sequence):
                left_process = self.process_map[left]
                right_process = self.process_map[right]
                plink(left_process, right_process)

    def _parse(self, filebuffer):
        raise NotImplementedError(
            'Inherit from this class and implement this method')


class InvestigationParser(AbstractParser):

    def __init__(self):
        self.isa = Investigation()
        self._investigation_prefix = 'Investigation'
        self._study_prefix = 'Study'
        self._term_accession_postfix = 'Term Accession Number'
        self._term_source_ref_postfix = 'Term Source REF'

    @staticmethod
    def _section_to_dict(section_rows):
        d = dict()
        for r in section_rows:
            label = next(iter(r), None)
            d[label] = r[1:]
        return d

    def _split_investigation_table(self, filebuffer):
        section_keywords = (
            'ONTOLOGY SOURCE REFERENCE', 'INVESTIGATION',
            'INVESTIGATION PUBLICATIONS', 'INVESTIGATION CONTACTS', 'STUDY',
            'STUDY DESIGN DESCRIPTORS', 'STUDY PUBLICATIONS', 'STUDY FACTORS',
            'STUDY ASSAYS', 'STUDY PROTOCOLS', 'STUDY CONTACTS')
        section_slices = []
        section_delimiters = []
        tabreader = csv.reader(filebuffer, delimiter='\t')
        for sec_index, row in enumerate(tabreader):
            label = next(iter(row), None)
            if label is not None and not label.strip().startswith('#'):
                if label in section_keywords:
                    section_delimiters.append(sec_index)
        filebuffer.seek(0)
        for this_sec_index, next_sec_index in self._pairwise(
                section_delimiters):
            section_slice = []
            sec_f = itertools.islice(filebuffer, this_sec_index, next_sec_index)
            secreader = csv.reader(sec_f, delimiter='\t')
            for row in secreader:
                section_slice.append(row)
            filebuffer.seek(0)
            section_slices.append(section_slice)
        sec_f = itertools.islice(filebuffer, section_delimiters[-1], None)
        section_slice = []
        secreader = csv.reader(sec_f, delimiter='\t')
        for row in secreader:
            section_slice.append(row)
        section_slices.append(section_slice)
        return section_slices

    @staticmethod
    def _parse_vertical_comments(section, section_obj):
        CommentGroup = namedtuple('CommentGroup', ['name', 'value'])
        regex_comment = re.compile('Comment\[(.*?)\]')
        comments = [(regex_comment.findall(x)[-1], section.get(x))
                    for x in (x for x in section.keys() if
                              regex_comment.match(x))]
        for comment in comments:
            comment_group = CommentGroup._make(comment)
            for value in getattr(comment_group, 'value'):
                section_obj.comments.append(
                    Comment(getattr(comment_group, 'name'), value))

    @staticmethod
    def _parse_horizontal_comments(section, section_obj_list):
        CommentGroup = namedtuple('CommentGroup', ['name', 'value'])
        regex_comment = re.compile('Comment\[(.*?)\]')
        comment_groups = [CommentGroup._make((regex_comment.findall(x)[-1],
                                             section.get(x)))
                          for x in (x for x in section.keys() if
                                    regex_comment.match(x))]
        for comment_group in comment_groups:
            for v, o in zip_longest(getattr(comment_group, 'value'),
                                    section_obj_list, fillvalue=''):
                o.comments.append(Comment(getattr(comment_group, 'name'), v))

    def _parse_ontology_source_reference_section(self, section):
        term_source_prefix = 'Term Source'
        names = section.get(' '.join([term_source_prefix, 'Name']), [])
        files = section.get(' '.join([term_source_prefix, 'File']), [])
        versions = section.get(' '.join([term_source_prefix, 'Version']), [])
        descriptions = section.get(
            ' '.join([term_source_prefix, 'Description']), [])
        for n, f, v, d in zip_longest(
                names, files, versions, descriptions, fillvalue=''):
            ontology_source = OntologySource(n, f, v, d)
            self.isa.ontology_source_references.append(ontology_source)
        self._parse_horizontal_comments(
            section, self.isa.ontology_source_references)

    def _parse_investigation_section(self, section, filename):
        identifier = section.get(
            ' '.join([self._investigation_prefix, 'Identifier']), [])
        title = section.get(
            ' '.join([self._investigation_prefix, 'Title']), [])
        description = section.get(
            ' '.join([self._investigation_prefix, 'Description']), [])
        submission_date = section.get(
            ' '.join([self._investigation_prefix, 'Submission Date']), [])
        public_release_date = section.get(
            ' '.join([self._investigation_prefix, 'Public Release Date']), [])
        i, t, d, s, p = next(zip_longest(
            identifier, title, description, submission_date,
            public_release_date, fillvalue=''), ('', '', '', '', ''))
        self.isa.identifier = i
        self.isa.title = t
        self.isa.description = d
        self.isa.submission_date = s
        self.isa.public_release_date = p
        self.isa.filename = filename
        self._parse_vertical_comments(section, self.isa)

    def _parse_study_section(self, section):
        identifier = section.get(
            ' '.join([self._study_prefix, 'Identifier']), [])
        title = section.get(
            ' '.join([self._study_prefix, 'Title']), [])
        description = section.get(
            ' '.join([self._study_prefix, 'Description']), [])
        submission_date = section.get(
            ' '.join([self._study_prefix, 'Submission Date']), [])
        public_release_date = section.get(
            ' '.join([self._study_prefix, 'Public Release Date']), [])
        filename = section.get(
            ' '.join([self._study_prefix, 'File Name']), [])
        i, t, d, s, p, f = next(zip_longest(
            identifier, title, description, submission_date,
            public_release_date, filename, fillvalue=''),
            ('', '', '', '', '', ''))
        study = Study('', f, i, t, d, s, p)
        self.isa.studies.append(study)
        self._parse_vertical_comments(section, study)

    def _parse_study_design_descriptors_section(self, section):
        design_type_prefix = ' '.join([self._study_prefix, 'Design Type'])
        design_type = section.get(design_type_prefix)
        design_type_term_accession = section.get(
            ' '.join([design_type_prefix, self._term_accession_postfix]), [])
        design_type_term_source_ref = section.get(
            ' '.join([design_type_prefix, self._term_source_ref_postfix]), [])
        for term, accession, source in zip_longest(
                design_type, design_type_term_accession,
                design_type_term_source_ref, fillvalue=''):
            if not all(x == '' for x in (term, accession, source)):
                annotation = OntologyAnnotation(term, source, accession)
                study = self.isa.studies[-1]
                study.design_descriptors.append(annotation)
        self._parse_horizontal_comments(
            section, self.isa.studies[-1].design_descriptors)

    def _parse_study_factors_section(self, section):
        factor_prefix = 'Factor'
        factor_type_prefix = ' '.join(
            [self._study_prefix, factor_prefix, 'Type'])
        factor_name = section.get(
            ' '.join([self._study_prefix, factor_prefix, 'Name']), [])
        factor_type = section.get(factor_type_prefix)
        factor_type_term_accession = section.get(
            ' '.join(
                [factor_type_prefix, self._term_accession_postfix]), [])
        factor_type_term_source_ref = section.get(
            ' '.join(
                [factor_type_prefix, self._term_source_ref_postfix]),
            [])
        for name, term, accession, source in zip_longest(
                factor_name, factor_type, factor_type_term_accession,
                factor_type_term_source_ref, fillvalue=''):
            if not all(x == '' for x in
                       (name, term, accession, source)):
                annotation = OntologyAnnotation(term, source, accession)
                factor = StudyFactor('', name, annotation)
                study = self.isa.studies[-1]
                study.factors.append(factor)
        self._parse_horizontal_comments(
            section, self.isa.studies[-1].factors)

    def _parse_study_assays_section(self, section):
        study_assay_prefix = ' '.join([self._study_prefix, 'Assay'])
        mt_prefix = ' '.join([study_assay_prefix, 'Measurement Type'])
        tt_prefix = ' '.join([study_assay_prefix, 'Technology Type'])
        measurement_type = section.get(mt_prefix)
        measurement_type_term_accession = section.get(
            ' '.join([mt_prefix, self._term_accession_postfix]), [])
        measurement_type_term_source_ref = section.get(
            ' '.join([mt_prefix, self._term_source_ref_postfix]), [])

        technology_type = section.get(tt_prefix, [])
        technology_type_term_accession = section.get(
            ' '.join([tt_prefix, self._term_accession_postfix]), [])
        technology_type_term_source_ref = section.get(
            ' '.join([tt_prefix, self._term_source_ref_postfix]), [])
        technology_platform = section.get(
            ' '.join([study_assay_prefix, 'Technology Platform']), [])
        assay_filename = section.get(
            ' '.join([study_assay_prefix, 'File Name']), [])
        for mt_term, mt_accession, mt_source, tt_term, tt_accession, \
            tt_source, p, f in zip_longest(
            measurement_type, measurement_type_term_accession,
            measurement_type_term_source_ref, technology_type,
            technology_type_term_accession,
            technology_type_term_source_ref, technology_platform,
            assay_filename, fillvalue=''):
            if not all(x == '' for x in (
                    mt_term, mt_accession, mt_source, tt_term, tt_accession,
                    tt_source, p, f)):
                mt = OntologyAnnotation(mt_term, mt_source, mt_accession)
                tt = OntologyAnnotation(tt_term, tt_source, tt_accession)
                assay = Assay(mt, tt, p, f)
                study = self.isa.studies[-1]
                study.assays.append(assay)
        self._parse_horizontal_comments(
            section, self.isa.studies[-1].assays)

    def _parse_study_protocols_section(self, section):
        protocol_prefix = ' '.join([self._study_prefix, 'Protocol'])
        protocol_type_prefix = ' '.join([protocol_prefix, 'Type'])
        parameters_name_prefix = ' '.join(
            [protocol_prefix, 'Parameters Name'])
        components_type_prefix = ' '.join(
            [protocol_prefix, 'Components Type'])
        protocol_name = section.get(
            ' '.join([protocol_prefix, 'Name']), [])
        protocol_type = section.get(protocol_type_prefix)
        protocol_type_accession = section.get(
            ' '.join([protocol_type_prefix, self._term_accession_postfix]),
            [])
        protocol_type_source_ref = section.get(
            ' '.join([protocol_type_prefix, self._term_source_ref_postfix]),
            [])
        protocol_description = section.get(
            ' '.join([protocol_prefix, 'Description']), [])
        protocol_uri = section.get(' '.join([protocol_prefix, 'URI']),
                                   [])
        protocol_version = section.get(
            ' '.join([protocol_prefix, 'Version']), [])
        parameters_names = section.get(parameters_name_prefix, [])
        parameters_names_term_accession = section.get(
            ' '.join(
                [parameters_name_prefix, self._term_accession_postfix]),
            [])
        parameters_names_term_source_ref = section.get(
            ' '.join([parameters_name_prefix,
                      self._term_source_ref_postfix]),
            [])
        components_names = section.get(
            ' '.join([protocol_prefix, 'Components Name']), [])
        components_types = section.get(components_type_prefix, [])
        components_types_term_accession = section.get(
            ' '.join(
                [components_type_prefix, self._term_accession_postfix]),
            [])
        components_types_term_source_ref = section.get(
            ' '.join([components_type_prefix,
                      self._term_source_ref_postfix]),
            [])
        for n, t, t_acc, t_src, d, u, v, pn, pn_acc, pn_src, cn, ct, \
            ct_acc, ct_src in zip_longest(
            protocol_name, protocol_type,
            protocol_type_accession,
            protocol_type_source_ref,
            protocol_description,
            protocol_uri, protocol_version,
            parameters_names,
            parameters_names_term_accession,
            parameters_names_term_source_ref,
            components_names,
            components_types,
            components_types_term_accession,
            components_types_term_source_ref, fillvalue=''):
            if not all(x == '' for x in (
                    n, t, t_acc, t_src, d, u, v, pn, pn_acc, pn_src, cn, ct,
                    ct_acc, ct_src)):
                t_ann = OntologyAnnotation(t, t_src, t_acc)
                protocol = Protocol('', n, t_ann, u, d, v)
                # parse Parameters
                for n, a, s in zip_longest(pn.split(';'),
                                           pn_acc.split(';'),
                                           pn_src.split(';'),
                                           fillvalue=''):
                    if not all(x == '' for x in (n, a, s)):
                        pn_ann = OntologyAnnotation(n, s, a)
                        parameter = ProtocolParameter('', pn_ann)
                        protocol.parameters.append(parameter)
                # parse Components
                for n, t, a, s in zip_longest(
                        cn.split(';'),
                        ct.split(';'),
                        ct_acc.split(';'),
                        ct_src.split(';'), fillvalue=''):
                    if not all(x == '' for x in (n, t, a, s)):
                        ct_ann = OntologyAnnotation(t, s, a)
                        component = ProtocolComponent('', n, ct_ann)
                        protocol.components.append(component)
                study = self.isa.studies[-1]
                study.protocols.append(protocol)
        self._parse_horizontal_comments(
            section, self.isa.studies[-1].protocols)

    def _parse_publications_section(self, section_label, section):
        investigation_or_study_prefix = \
            'Investigation' if 'INVESTIGATION' in next(
                iter(section_label)) else 'Study'
        publication_prefix = 'Publication'
        status_prefix = 'Status'
        pubmed_id = section.get(' '.join(
            [investigation_or_study_prefix, 'PubMed ID']), [])
        doi = section.get(' '.join(
            [investigation_or_study_prefix, publication_prefix,
             'DOI']), [])
        author_list = section.get(' '.join(
            [investigation_or_study_prefix, publication_prefix,
             'Author List']), [])
        title = section.get(' '.join(
            [investigation_or_study_prefix, publication_prefix,
             'Title']), [])
        status = section.get(' '.join(
            [investigation_or_study_prefix, publication_prefix,
             status_prefix]))
        status_accession = section.get(' '.join(
            [investigation_or_study_prefix, publication_prefix,
             status_prefix, self._term_accession_postfix]))
        status_term_source = section.get(' '.join(
            [investigation_or_study_prefix, publication_prefix,
             status_prefix, self._term_source_ref_postfix]))
        for p, d, a, t, s, s_acc, s_src in zip_longest(
                pubmed_id, doi, author_list, title, status,
                status_accession, status_term_source, fillvalue=''):
            if not all(x == '' for x in (p, d, a, t, s, s_acc, s_src)):
                annotation = OntologyAnnotation(s, s_src, s_acc)
                publication = Publication(p, d, a, t, annotation)
                if next(iter(self.isa.studies), None) is None:
                    self.isa.publications.append(publication)
                else:
                    study = self.isa.studies[-1]
                    study.publications.append(publication)
        if next(iter(self.isa.studies), None) is None:
            self._parse_horizontal_comments(
                section, self.isa.publications)
        else:
            self._parse_horizontal_comments(
                section, self.isa.studies[-1].publications)

    def _parse_contacts_section(self, section_label, section):
        investigation_or_study_prefix = \
            'Investigation' if 'INVESTIGATION' in next(
                iter(section_label)) else 'Study'
        person_prefix = 'Person'
        roles_prefix = 'Roles'
        last_name = section.get(' '.join(
            [investigation_or_study_prefix, person_prefix,
             'Last Name']), [])
        first_name = section.get(' '.join(
            [investigation_or_study_prefix, person_prefix,
             'First Name']), [])
        mid_initials = section.get(' '.join(
            [investigation_or_study_prefix, person_prefix,
             'Mid Initials']), [])
        email = section.get(' '.join(
            [investigation_or_study_prefix, person_prefix, 'Email']),
            [])
        phone = section.get(' '.join(
            [investigation_or_study_prefix, person_prefix, 'Phone']),
            [])
        fax = section.get(' '.join(
            [investigation_or_study_prefix, person_prefix, 'Fax']), [])
        address = section.get(' '.join(
            [investigation_or_study_prefix, person_prefix, 'Address']),
            [])
        affiliation = section.get(' '.join(
            [investigation_or_study_prefix, person_prefix,
             'Affiliation']), [])

        roles = section.get(' '.join(
            [investigation_or_study_prefix, person_prefix,
             roles_prefix]), [])
        roles_term_accession = section.get(' '.join(
            [investigation_or_study_prefix, person_prefix,
             roles_prefix, self._term_accession_postfix]), [])
        roles_term_source_ref = section.get(' '.join(
            [investigation_or_study_prefix, person_prefix,
             roles_prefix, self._term_source_ref_postfix]), [])

        for l, fn, m, e, p, fx, ad, af, rs, r_accs, r_srcs in \
                zip_longest(last_name, first_name,
                            mid_initials, email, phone, fax,
                            address, affiliation, roles,
                            roles_term_accession,
                            roles_term_source_ref,
                            fillvalue=''):
            if not all(x == '' for x in (
                    l, fn, m, e, p, fx, ad, af, rs, r_accs, r_srcs)):
                roles_list = []
                for r, r_acc, r_src in zip_longest(
                        rs.split(';'), r_accs.split(';'),
                        r_srcs.split(';'), fillvalue=''):
                    if not all(x == '' for x in (r, r_acc, r_src)):
                        r_ann = OntologyAnnotation(r, r_src, r_acc)
                        roles_list.append(r_ann)
                person = Person(l, fn, m, e, p, fx, ad, af, roles_list)
                if next(iter(self.isa.studies), None) is None:
                    self.isa.contacts.append(person)
                else:
                    study = self.isa.studies[-1]
                    study.contacts.append(person)
        if next(iter(self.isa.studies), None) is None:
            self._parse_horizontal_comments(
                section, self.isa.contacts)
        else:
            self._parse_horizontal_comments(
                section, self.isa.studies[-1].contacts)

    def _parse(self, filebuffer):
        section_slices = self._split_investigation_table(filebuffer)
        for section_slice in section_slices:
            section_label = next(iter(section_slice))
            section = self._section_to_dict(section_slice)
            if 'ONTOLOGY SOURCE REFERENCE' in section_label:
                self._parse_ontology_source_reference_section(section)
            if 'INVESTIGATION' in section_label:
                self._parse_investigation_section(section, os.path.basename(
                    getattr(filebuffer, 'name', '')))
            if 'STUDY' in section_label:
                self._parse_study_section(section)
            if 'STUDY DESIGN DESCRIPTORS' in section_label:
                self._parse_study_design_descriptors_section(section)
            if 'STUDY FACTORS' in section_label:
                self._parse_study_factors_section(section)
            if 'STUDY ASSAYS' in section_label:
                self._parse_study_assays_section(section)
            if 'STUDY PROTOCOLS' in section_label:
                self._parse_study_protocols_section(section)
            if any(x in section_label for x in (
                    'INVESTIGATION PUBLICATIONS', 'STUDY PUBLICATIONS')):
                self._parse_publications_section(section_label, section)
            if any(x in section_label for x in (
                    'INVESTIGATION CONTACTS', 'STUDY CONTACTS')):
                self._parse_contacts_section(section_label, section)
