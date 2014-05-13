"""Module to implement ORM to the ome database"""

from types import MethodType
from os import system

from sqlalchemy.orm import sessionmaker, relationship, aliased
from sqlalchemy.orm.session import Session as _SA_Session
from sqlalchemy import Table, MetaData, create_engine,Column, Integer, \
    String, Float, ForeignKey, and_, or_, not_, distinct, select
from sqlalchemy.schema import UniqueConstraint
from sqlalchemy.ext.declarative import declarative_base
import om_dev.lib.settings as settings
import pymongo


engine = create_engine("postgresql://%s:%s@%s/%s" %
    (settings.user, settings.password, settings.host, settings.dev_database))
Base = declarative_base(bind=engine)
metadata = MetaData(bind=engine, schema=settings.schema)

connection = pymongo.Connection()
#omics_database = connection.omics_database
omics_database = connection.omics_database3

class id2otherid(Base):
    __tablename__ = "id2otherid"
    
    id = Column(Integer, primary_key=True)
    other_id = Column(String(100), primary_key=True)
    type = Column(String(25))
    data_source_id = Column(Integer, ForeignKey('data_source.id'))
    data_source = relationship("DataSource")
    
    __table_args__ = (UniqueConstraint('id','other_id'),{})

    def __repr__(self):
        return "%s in (%s)" % (self.other_id, str(self.data_source))
    
    def __init__(self, id, other_id, type, data_source_id):
        self.id = id
        self.other_id = other_id
        self.type = type
        self.data_source_id = data_source_id
        

class DataSource(Base):
    __tablename__ = 'data_source'
    
    id = Column(Integer, primary_key=True)
    name = Column(String(100))
    lab = Column(String(100))
    institution = Column(String(100))
    #data_sets = relationship("DataSet")

    __table_args__ = (UniqueConstraint('name'),{})
    
    def __repr__(self):
        return "Data Source %s (#%d)" % (self.name, self.id)
    
    def __repr__dict__(self):
        return {"name":self.name,"wid":self.id,"values":{"lab":self.lab,"institution":self.institution}}
    
    def __repr__json__(self):
        return json.dumps(self.__repr__dict__())
    
    def __init__(self, name, lab=None, institution=None):
        self.name = name
        self.lab = lab
        self.institution = institution
        

class GenomeRegion(Base):
    __tablename__ = 'genome_region'    
    
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(10))
    leftpos = Column(Integer, nullable=False)
    rightpos = Column(Integer, nullable=False)
    strand = Column(String(1), nullable=False)
    type = Column(String(20))
    
    __table_args__ = (UniqueConstraint('leftpos','rightpos','strand'),{})

    __mapper_args__ = {'polymorphic_identity': 'genome_region', 
                       'polymorphic_on': type
                      }
    
    def __repr__(self):
        return "GenomeRegion: %d-%d (%s)" % \
                (self.leftpos, self.rightpos, self.strand)
                
    def __init__(self, leftpos, rightpos, strand, name=None):
        self.leftpos = leftpos
        self.rightpos = rightpos
        self.strand = strand
        self.name = name
        
        
def make_table(table_name):
    """function to create a table with the default parameters"""
    return Table(table_name, metadata, autoload=True)


class _Session(_SA_Session):
    """an sqlalchemy session object to interact with the OME database

    This object can used to make queries against the ome database. For
    example, a query without using any ORM looks like this
    >>> session = Session()
    >>> session.execute("SELECT name from genes where bnum='b0001'").fetchone()
    (u'thrL',)
    Using the sqlalchemy ORM gives more descriptive objects. For example:
    >>> b0001 = session.query(Gene).filter(Gene.bnum=="b0001").first()
    >>> b0001.name
    u'thrL'
    Raw queries which return ORM objects are also possible:
    >>> sql_statement = "SELECT * from genes where bnum='b0001'"
    >>> b0001 = session.query(Gene).from_statement(sql_statement).first()
    >>> b0001.name
    u'thrL'

    The Session will automatically set the search_path to settings.schema
    """
    
    
    def __init__(self, *args, **kwargs):
        super(_Session, self).__init__(*args, **kwargs)
        self.execute("set search_path to %s;" % (settings.schema))
        self.commit()
        self.get_or_create = MethodType(get_or_create, self)
        #self.search_by_synonym = MethodType(search_by_synonym, self)


    def __repr__(self):
        return "OME session %d" % (self.__hash__())


def get_or_create(session, class_type, **kwargs):
    """gets an object using filter_by on the unique kwargs. If no such object
    is found in the database, a new one will be created which satisfies
    these constraints. This is why every class that wants to use this
    method to be instantiated needs to have a UniqueConstraint defined."""
    
    for constraint in list(class_type.__table_args__):
        if constraint.__class__.__name__ == 'UniqueConstraint':
            unique_cols = constraint.columns.keys()
    #print class_type.__name__
    #print unique_cols
    #print kwargs
    try:
        result = session.query(class_type).filter_by(**{k: kwargs[k] for k in unique_cols}).one()
    except:
        session.add(class_type(**kwargs))
        session.commit()
        result = session.query(class_type).filter_by(**{k: kwargs[k] for k in unique_cols}).one()
    return result


def update(session, object, **kwargs):
    """Ideally this would only search on the primary key columns so
    that an update could be made in one call. However, its not currently
    clear how to do that so necessary to pass in the actual object and
    update following a call to get_or_create() There is probably some
    way to do this with class_mapper but its hard right now
    """
    #result = session.query(class_type).filter_by(**kwargs).first()
    #result = session.query(class_type).filter_by(name=kwargs['name']).first()
    #if result is None: return
    
    for key,value in kwargs.iteritems(): 
        setattr(object,key,value) 
    session.add(object)
    session.commit()
    
    return object
        
        
Session = sessionmaker(bind=engine, class_=_Session)


if __name__ == "__main__":
    session = Session()
    
    