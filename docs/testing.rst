.. _testing:

Testing
=======

Each module within the source code has its own module for unit testing.
Please see :ref:`logging_tests`,
:ref:`data_handling_tests`,
and :ref:`multithreading_tests` for more information.
Additionally,
there is a class-based "fixture" - 
well it's not really a fixture,
pytest doesn't allow class based fixtures - 
for creating custom temporary file.
Please see :ref:`custom_temp_file`.
There are long term plans to develop integration tests for the whole pipeline.

.. toctree::
   :hidden:
   :maxdepth: 3

   logging_tests
   data_handling_tests
   multithreading_tests
   custom_temp_file
